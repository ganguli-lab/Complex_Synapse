function [ beta ] = BWbetaN( eta,modelobj,simobj,updater )
%beta=BWBETA(eta,modelobj,simobj) normalised backward variables for Baum-Welch algorithm
%   beta     = normalised backward variables
%   eta      = normalisation factor
%   modelobj = SynapseIdModel
%   simobj   = SynapsePlastSeq

%   t        = time-step after which we want forward variables

error(CheckSize(modelobj,@isscalar));
error(CheckSize(simobj,@isscalar));
% error(CheckSize(t,@isscalar));
% error(CheckValue(t,@isint));

% M=cat(3,modelobj.M{simobj.potdep});
% outProj=cat(3,modelobj.outProj{simobj.readouts});

beta=ones(modelobj.NumStates,simobj.NumT);

if nargin>=5
    for i=simobj.NumT:-1:2
        beta(:,i-1) = squeeze(updater(:,:,i-1)) * beta(:,i) * eta(i);
    end
else
    for i=simobj.NumT:-1:2
        beta(:,i-1) = modelobj.M{simobj.potdep(i-1)} * modelobj.outProj{simobj.readouts(i)} * beta(:,i) * eta(i);
    end
end

end

