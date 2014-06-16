function [ alpha,eta,updater ] = BWalphaN( modelobj,simobj )
%[alpha,eta]=BWALPHAN(modelobj,simobj) normalised forward variables for Baum-Welch algorithm
%   alpha    = normalised forward variables
%   eta      = normalisation factor
%   modelobj = SynapseIdModel
%   simobj   = SynapsePlastSeq

%   t        = time-step up to which we want forward variables

error(CheckSize(modelobj,@isscalar));
error(CheckSize(simobj,@isscalar));
% error(CheckSize(t,@isscalar));
% error(CheckValue(t,@isint));


alpha=zeros(simobj.NumT,modelobj.NumStates);
eta=zeros(simobj.NumT,1);
alpha(1,:) = modelobj.Initial * modelobj.outProj{simobj.readouts(1)};
eta(1) = 1 / sum(alpha(1,:));
alpha(1,:) = alpha(1,:) * eta(1);

if nargout >= 3
    M=cat(3,modelobj.M{simobj.potdep(1:end-1)});
    outProj=cat(3,modelobj.outProj{simobj.readouts(2:end)});
    updater=M;
    for i=2:simobj.NumT
        updater(:,:,i-1) = squeeze(M(:,:,i-1)) * squeeze(outProj(:,:,i-1));
        alpha(i,:) = alpha(i-1,:) * squeeze(updater(:,:,i-1));
        eta(i) = 1 / sum(alpha(i,:));
        alpha(i,:) = alpha(i,:) * eta(i);
    end
else
    for i=2:simobj.NumT
        alpha(i,:) = alpha(i-1,:) * modelobj.M{simobj.potdep(i-1)} * modelobj.outProj{simobj.readouts(i)};
        eta(i) = 1 / sum(alpha(i,:));
        alpha(i,:) = alpha(i,:) * eta(i);
    end
end%if nargout

end

