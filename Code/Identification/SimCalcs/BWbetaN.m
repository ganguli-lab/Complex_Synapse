function [ beta ] = BWbetaN( eta,t,modelobj,simobj )
%beta=BWBETA(eta,t,modelobj,simobj) normalised backward variables for Baum-Welch algorithm
%   beta     = normalised backward variables
%   eta      = normalisation factor
%   t        = time-step after which we want forward variables
%   modelobj = SynapseIdModel
%   simobj   = SynapsePlastSeq


error(CheckSize(t,@isscalar));
error(CheckValue(t,@isint));

beta=ones(length(modelobj.Initial),length(simobj.readouts)-t+1);

for i=length(simobj.readouts):-1:t+1
    beta(:,i-t) = modelobj.M{simobj.potdep(i-1)} * modelobj.outProj{simobj.readouts(i)} * beta(:,i-t+1) * eta(i);
end

end

