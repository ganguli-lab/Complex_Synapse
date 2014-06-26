function [ optimValues,updaterfn,extraArgs ] = FitSynapseUtility_setup( options,guessmodel,simobj,init )
%[optimValues,updaterfn,extraArgs]=FITSYNAPSEUTILITY_SETUP(options,guessmodel,simobj,init)
%initial setup of variables for optimiser
%   optimValues = struct with information about the current state of the optimiser
%   simobj      = vector of SynapsePlastSeq with data to be fit
%   updaterfn   = handle of function that'll perform the updates
%   extraArgs   = cell array of extra arguments to pass to updaterfn
%   options     = struct of options (see SynapseOptimset)
%   guessmodel  = SynapseIdModel with initial guess of the model
%   init        = 'Init' if we're only updating fitmodel.Initial (optional)


existsAndDefault('init','');

%information about the current state of the optimiser
optimValues=struct('iteration',0,'funcCount',0,...
    'procedure',[options.Algorithm ',' options.Weighter ',' options.Penaliser],...
    'NumStates',guessmodel.NumStates,'NumPlast',guessmodel.NumPlast,...
    'model',[],'fval',[],'prev',[],'truth',[],'stateProb',{{}});

%the SynapseIdModel where we'll store the fit
optimValues.model=guessmodel.Sort;
optimValues.fval=HMMloglike(optimValues.model,simobj)+SynapsePrior(optimValues.model,options);

%information about the previous state and current state of the optimiser relative to the previous state
optimValues.prev=struct('model',[],'fval',[],'dfval',[],'dInitial',[],'dM',[]);
%information about the current state of the optimiser relative to ground truth
if ~isempty(options.GroundTruth)
    optimValues.truth=optimValues.prev;
    optimValues.truth.model=options.GroundTruth;
    optimValues.truth.fval=HMMloglike(options.GroundTruth,simobj)+SynapsePrior(options.GroundTruth,options);
    optimValues.truth.dfval=optimValues.truth.fval-optimValues.fval;
end

updaterfn=str2func([options.Penaliser 'penalty']);
extraArgs=[{'Algorithm',[init options.Algorithm],'Weighter',options.Weighter,'Penalty',options.Penalty},options.ExtraParams];



end

