function [ beta ] = BWbetaN( eta_up,modelobj,simobj )
%beta=BWBETA(eta,modelobj,simobj) normalised backward variables for Baum-Welch algorithm
%   beta     = normalised backward variables
%   eta_up   = normalisation factor or updater matrices
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

if ismatrix(eta_up)
    for i=simobj.NumT:-1:2
        beta(:,i-1) = modelobj.M{simobj.potdep(i-1)} * modelobj.outProj{simobj.readouts(i)} * beta(:,i) * eta_up(i);
    end
else
    siz=modelobj.NumStates*[1 1];
    for i=simobj.NumT:-1:2
        beta(:,i-1) = reshape(eta_up(:,:,i-1),siz) * beta(:,i);
    end
end

end

