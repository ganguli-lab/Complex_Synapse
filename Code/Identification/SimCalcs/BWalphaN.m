function [ alpha,eta ] = BWalphaN( t,modelobj,simobj )
%[alpha,eta]=BWALPHAN(t,modelobj,simobj) normalised forward variables for Baum-Welch algorithm
%   alpha    = normalised forward variables
%   eta      = normalisation factor
%   t        = time-step up to which we want forward variables
%   modelobj = SynapseIdModel
%   simobj   = SynapsePlastSeq


error(CheckSize(t,@isscalar));
error(CheckValue(t,@isint));

alpha=zeros(t,length(modelobj.Initial));
eta=zeros(t,1);
alpha(1,:) = modelobj.Initial * modelobj.outProj{simobj.readouts(1)};
eta(1) = 1 / sum(alpha(1,:));
alpha(1,:) = alpha(1,:) * eta(1);

for i=2:t
    alpha(i,:) = alpha(i-1,:) * modelobj.M{simobj.potdep(i-1)} * modelobj.outProj{simobj.readouts(i)};
    eta(i) = 1 / sum(alpha(i,:));
    alpha(i,:) = alpha(i,:) * eta(i);
end

end

