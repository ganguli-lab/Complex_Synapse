function [ optimValues,exitflag,msg,ME ] = FitSynapseUtility_update( optimValues,options,simobj,updaterfn,extraArgs )
%[optimValues,exitflag,msg]=FITSYNAPSEUTILITY_UPDATE(optimValues,simobj,updaterfn,extraArgs)
%perform update of currrent model fit
%   optimValues = struct with information about the current state of the optimiser
%   options     = struct of options (see SynapseOptimset)
%   simobj      = vector of SynapsePlastSeq with data to be fit
%   updaterfn   = handle of function that'll perform the updates
%   extraArgs   = cell array of extra arguments to pass to updaterfn
%   exitflag    = integer describing reason for exiting
%   msg         = string describing reason for exiting
%   ME          = MException

exitflag=0;
msg=['Exceeded max iterations: ' int2str(options.MaxIter)];
ME=[];

%store previous model
    %
    optimValues.iteration=optimValues.iteration+1;
    optimValues.prev.model=optimValues.model;
    optimValues.prev.fval=optimValues.fval;
    %
%perform the update
    %
    try
    [optimValues.model,optimValues.fval,optimValues.stateProb]=updaterfn(optimValues.prev.model,simobj,...
        extraArgs{:});
    optimValues.model=optimValues.model.Sort(options.fp);
    catch ME
        exitflag=-1;
        msg=['Error: ' ME.message];
    end


end

