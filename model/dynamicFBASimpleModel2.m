function out = dynamicFBASimpleModel2()
% This script generates an FBA model for the simple network used in the
% INSA course, and allows the model to be analyzed by means of the classical FBA
% method, dynamic FBA, or dynamic regulatory FBA (the latter using an ugly
% hack).
% In order to run the script, the COBRA Matlab toolbox needs to be
% installed
%
% HdJ 16/3/17

clear all;
close all;

% Initialize COBRA v2 toolbox
initCobraToolbox;
solverOK = changeCobraSolver('glpk','LP')

% Define run options
doFBA = 0;
doDynamicFBA = 1;
doRegulatoryDynamicFBA = 0;


% Define model
model = createModel({},{},{},[],[],[]);
model = addReaction(model,'EX_S1','S1[e] <=> '); model = changeRxnBounds(model,'EX_S1',-1000,'l'); model = changeRxnBounds(model,'EX_S1',1000,'u');
model = addReaction(model,'Import_S1','S1[e] <=> X1'); model = changeRxnBounds(model,'Import_S1',-10,'l'); model = changeRxnBounds(model,'Import_S1',10,'u');
model = addReaction(model,'EX_S2','S2[e] <=> '); model = changeRxnBounds(model,'EX_S2',-1000,'l'); model = changeRxnBounds(model,'EX_S2',1000,'u');
model = addReaction(model,'Import_S2','S2[e] <=> X2'); model = changeRxnBounds(model,'Import_S2',-10,'l'); model = changeRxnBounds(model,'Import_S2',10,'u');
model = addReaction(model,'MetReaction_1','X1 <=> 4 M'); model = changeRxnBounds(model,'MetReaction_1',-1000,'l'); model = changeRxnBounds(model,'MetReaction_1',1000,'u');
model = addReaction(model,'MetReaction_2','X2 <=> M'); model = changeRxnBounds(model,'MetReaction_2',-1000,'l'); model = changeRxnBounds(model,'MetReaction_2',1000,'u');
model = addReaction(model,'Biomass','10 M ->'); model = changeRxnBounds(model,'Biomass',0,'l'); model = changeRxnBounds(model,'Biomass',1000,'u'); 

% Initialize objective function
model.c = [0 0 0 0 0 0 1]';

% Define dynamic FBA
substrateRxns = {'EX_S1', 'EX_S2'};
initConcentrations = [50 20]; initBiomass = 0.1;
timeStep = .02; nSteps = 150;
plotRxns = {'EX_S1','EX_S2','Biomass'};

% Run classical FBA
if doFBA
    optimizeCbModel(model,'max')
end

% Run dynamic FBA
if doDynamicFBA
    dynamicFBA(model,substrateRxns,initConcentrations,initBiomass,timeStep,nSteps,plotRxns);
end

% Run regulatory dynamic FBA, using ugly hack in dynamicFBA function in
% COBRA
if doRegulatoryDynamicFBA
    dynamicFBA_regulation(model,substrateRxns,initConcentrations,initBiomass,timeStep,nSteps,plotRxns);
end


end