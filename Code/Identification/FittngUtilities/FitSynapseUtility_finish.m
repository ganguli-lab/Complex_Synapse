function [ fval,output ] = FitSynapseUtility_finish( optimValues,options,msg,ME )
%[fval,output]=FitSynapseUtility_finish(optimValues,options,msg)
%create ouputs after optimiser terminates
%   fval   = log likelihood of fitmodel
%   output = struct that contains information about the optimization 
%   optimValues = struct with information about the current state of the optimiser
%   options     = struct of options (see SynapseOptimset)
%   msg         = string describing reason for exiting
%   ME          = MException

fval=optimValues.fval;
output=struct('message',msg,'ME',ME,...
    'algortihm',options.Algorithm,'weighter',options.Weighter,'penaliser',options.Penaliser,...
    'iterations',optimValues.iteration,'prev',optimValues.prev,'truth',optimValues.truth);

end

