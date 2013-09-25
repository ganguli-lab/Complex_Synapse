function [ loglike ] = HMMloglike( modelobj,simobj )
%like=HMMLOGLIKE(modelobj,simobj) log likelihood of outputs for Hidden-Markov-Model
%   loglike  = log likelihood
%   modelobj = SynapseIdModel
%   simobj   = SynapsePlastSeq



[~,eta]=BWalphaN(length(simobj.readouts),modelobj,simobj);
loglike = -sum(log(eta));

end

