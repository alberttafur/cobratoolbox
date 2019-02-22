function [dFBAsol] = dynamicFBA(model, substrateRxns, initConcentrations, plotRxns, varargin)
% Performs dynamic FBA simulation using the static optimization approach
%
% USAGE:
%
%    [dFBAsol] = dynamicFBA(model, substrateRxns, initConcentrations, initBiomass, timeStep, nSteps, plotRxns, exclUptakeRxns)
%
% INPUTS:
%    model:                  COBRA model structure
%    substrateRxns:          List of exchange reaction names for substrates
%                            initially in the media that may change (e.g. not
%                            h2o or co2)
%    initConcentrations:     Initial concentrations of substrates (in the same
%                            structure as `substrateRxns`)
%    plotRxns:               Reactions to be plotted (Default = {'EX_glc(e)', 'EX_ac(e)', 'EX_for(e)'})
%
% OPTIONAL INPUTS:
%    initBiomass:            Initial biomass (must be non zero)
%                            Default = 0.033;
%    timeStep:               Time step size
%                            Default = 0.1;
%    nSteps:                 Maximum number of time steps
%                            Default = 500;
%    exclUptakeRxns:         List of uptake reactions whose substrate concentrations do not change
%                            (Default = {'EX_co2(e)', 'EX_o2(e)', 'EX_h2o(e)', 'EX_h(e)'})
%    savePlot:               (double or boolean) boolean for saving plot in a file
%                            Default = false
%    showPlot:               (double or boolean) boolean for showing the plot
%                            Default = false
%    fileName:               (char) name of the file where the plot is saved.
%                            Default = product
%    outputFolder:           (char) name of the folder where files are saved
%                            Default = 'Results'
%
% OUTPUTS:
%    dFBAsol:    Matrix of time points, biomass values, and extracellular metabolite concentrations
%
% If no initial concentration is given for a substrate that has an open
% uptake in the model (i.e. `model.lb < 0`) the concentration is assumed to
% be high enough to not be limiting. If the uptake rate for a nutrient is
% calculated to exceed the maximum uptake rate for that nutrient specified
% in the model and the max uptake rate specified is > 0, the maximum uptake
% rate specified in the model is used instead of the calculated uptake
% rate.
%
% NOTE:
%
%    The dynamic FBA method implemented in this function is essentially
%    the same as the method described in
%    [`Varma, A., and B. O. Palsson. Appl. Environ. Microbiol. 60:3724 (1994)`].
%    This function does not implement the dynamic FBA using dynamic optimization approach
%    described in [`Mahadevan, R. et al. Biophys J, 83:1331-1340 (2003)`].
%
% .. Author: - Markus Herrgard 8/22/06
%            - Bug fix to plot, to save and export data, Albert Tafur Rangel, 2/22/19, Grupo de Diseño de Productos y Procesos, Universidad de Los Andes, Colombia.

parser = inputParser();
parser.addRequired('model', @(x) isstruct(x) && isfield(x, 'S') && isfield(model, 'rxns')...
    && isfield(model, 'mets') && isfield(model, 'lb') && isfield(model, 'ub') && isfield(model, 'b')...
    && isfield(model, 'c'));
parser.addRequired('substrateRxns', @(x)  iscell(x) && ~isempty(x));
parser.addRequired('initConcentrations', @(x) isnumeric(x) && ~isempty(x));
parser.addRequired('plotRxns', @(x)  iscell(x) && ~isempty(x));
parser.addParamValue('initBiomass', 0.033, @isnumeric);
parser.addParamValue('timeStep', 0.1, @isnumeric);
parser.addParamValue('nSteps', 500, @isnumeric);
parser.addParamValue('exclUptakeRxns', {'EX_co2(e)', 'EX_o2(e)', 'EX_h2o(e)', 'EX_h(e)'}, @(x)  iscell(x) && ~isempty(x));
parser.addParamValue('savePlot', 0, @(x) isnumeric(x) || islogical(x));
parser.addParamValue('showPlot', 0, @(x) isnumeric(x) || islogical(x));
parser.addParamValue('fileName', 'dFBA', @(x) ischar(x))
parser.addParamValue('outputFolder', '', @(x) ischar(x))

parser.parse(model, substrateRxns, initConcentrations, plotRxns, varargin{:})
model = parser.Results.model;
substrateRxns = parser.Results.substrateRxns;
initConcentrations = parser.Results.initConcentrations;
plotRxns = parser.Results.plotRxns;
initBiomass = parser.Results.initBiomass;
timeStep = parser.Results.timeStep;
nSteps = parser.Results.nSteps;
exclUptakeRxns = parser.Results.exclUptakeRxns;
savePlot = parser.Results.savePlot;
showPlot = parser.Results.showPlot;
fileName = parser.Results.fileName;
outputFolder = parser.Results.outputFolder;


global WAITBAR_TYPE

% Uptake reactions whose substrate concentrations do not change
if (nargin < 8)
    exclUptakeRxns = {'EX_co2(e)','EX_o2(e)','EX_h2o(e)','EX_h(e)'};
end

% Find exchange rxns
excInd = findExcRxns(model,false);
excInd = excInd & ~ismember(model.rxns,exclUptakeRxns);
excRxnNames = model.rxns(excInd);
length(excRxnNames)
% Figure out if substrate reactions are correct
missingInd = find(~ismember(substrateRxns,excRxnNames));
if (~isempty(missingInd))
    for i = 1:length(missingInd)
        fprintf('%s\n',substrateRxns{missingInd(i)});
    end
    error('Invalid substrate uptake reaction!');
end

% Initialize concentrations
[~, substrateMatchInd] = ismember(substrateRxns,excRxnNames);
concentrations = zeros(length(excRxnNames),1);
concentrations(substrateMatchInd) = initConcentrations;

% Deal with reactions for which there are no initial concentrations
originalBound = -model.lb(excInd);
noInitConcentration = (concentrations == 0 & originalBound > 0);
concentrations(noInitConcentration) = 1000;

biomass = initBiomass;

% Initialize bounds
uptakeBound =  concentrations/(biomass*timeStep);

% Make sure bounds are not higher than what are specified in the model
aboveOriginal = (uptakeBound > originalBound) & (originalBound > 0);
uptakeBound(aboveOriginal) = originalBound(aboveOriginal);
model.lb(excInd) = -uptakeBound;

concentrationMatrix = sparse(concentrations);
biomassVec = biomass;
timeVec(1) = 0;

fprintf('Step number\tBiomass\n');
showprogress(0,'Dynamic FBA analysis in progress ...');
for stepNo = 1:nSteps
    % Run FBA
    sol = optimizeCbModel(model,'max','one');
    mu = sol.f;
    if (sol.stat ~= 1 || mu == 0)
        fprintf('\nNo feasible solution - nutrients exhausted. Biomass:\t %f\n', biomass);
        break;
    end
    uptakeFlux = sol.x(excInd);
    biomass = biomass*exp(mu*timeStep);
    %biomass = biomass*(1+mu*timeStep);
    biomassVec(end+1) = biomass;

    % Update concentrations
    concentrations = concentrations - uptakeFlux/mu*biomass*(1-exp(mu*timeStep));
    %concentrations = concentrations + uptakeFlux*biomass*timeStep;
    concentrations(concentrations <= 0) = 0;
    concentrationMatrix(:,end+1) = sparse(concentrations);

    % Update bounds for uptake reactions
    uptakeBound =  concentrations/(biomass*timeStep);
    % This is to avoid any numerical issues
    uptakeBound(uptakeBound > 1000) = 1000;
    % Figure out if the computed bounds were above the original bounds
    aboveOriginal = (uptakeBound > originalBound) & (originalBound > 0);
    % Revert to original bounds if the rate was too high
    uptakeBound(aboveOriginal) = originalBound(aboveOriginal);
    uptakeBound(abs(uptakeBound) < 1e-9) = 0;

    model.lb(excInd) = -uptakeBound;

    if WAITBAR_TYPE ~= 1
        fprintf('%d\t%f\n',stepNo,biomass);
    end
    showprogress(stepNo/nSteps);
    timeVec(stepNo+1) = stepNo*timeStep;
end

selNonZero = any(concentrationMatrix>0,2);
concentrationMatrix = concentrationMatrix(selNonZero,:);
excRxnNames = excRxnNames(selNonZero);
selPlot = ismember(excRxnNames,plotRxns);

%create a table with the matrix concentration for each EX metbolite and Biomass
dFBAsol = [table(timeVec.', biomassVec.','VariableNames',{'Time','Biomass'}),...
    array2table(concentrationMatrix.', 'VariableNames',excRxnNames.')];

f = figure;
if ~showPlot
    set(gcf, 'Visible', 'Off');
end
%plot dFBA solution for the reactions selected
yyaxis left 
plot(timeVec,biomassVec)
ylabel('gDW / L')
xlabel('Time')
yyaxis right 
plot(timeVec,concentrationMatrix(selPlot,:))
ylabel('mmol / L')
legend(vertcat('biomass', strrep(excRxnNames(selPlot),'EX_',''))) 

if savePlot
    %directory change
    currentDirectory = pwd;
    if ~isempty(outputFolder)
        cd(outputFolder);
    else
        fullPath = which('tutorial_FBA');
        folder = fileparts(fullPath);
        mkdir([folder filesep 'FBAresults']);
        cd([folder filesep 'FBAresults']);
    end

    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 20 10]);
    saveas(gcf,[fileName '.png'])
    if ~showPlot
        close(f);
    end
    cd(currentDirectory);
end






