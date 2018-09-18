function [ loglike ] = HMMloglike( modelobj,simobj )
%like=HMMLOGLIKE(modelobj,simobj) log likelihood of outputs for Hidden-Markov-Model
%   loglike  = log likelihood
%   modelobj = SynapseIdModel
%   simobj   = SynapsePlastSeq


% error(CheckSize(modelobj,@isscalar));

if isscalar(simobj)
    [~,eta]=BWalphaN(modelobj,simobj);
    loglike = -sum(log(eta));
else
    loglike=0;
    for i=1:numel(simobj)
        loglike=loglike+HMMloglike(modelobj,simobj(i));
    end
end

end

