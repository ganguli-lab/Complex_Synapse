function [ alpha,eta,updater ] = BWalphaN( modelobj,simobj )
%[alpha,eta,updater]=BWALPHAN(modelobj,simobj) normalised forward variables for Baum-Welch algorithm
%   alpha    = normalised forward variables
%   eta      = normalisation factor
%   updater  = updater matrices (M*outProj*eta)
%   modelobj = SynapseIdModel
%   simobj   = SynapsePlastSeq


% error(CheckSize(modelobj,@isscalar));
% error(CheckSize(simobj,@isscalar));


alpha=zeros(simobj.NumT,modelobj.NumStates);
eta=zeros(simobj.NumT,1);

alpha(1,:) = modelobj.Initial * modelobj.outProj{simobj.readouts(1)};
eta(1) = 1 / sum(alpha(1,:));
alpha(1,:) = alpha(1,:) * eta(1);

M=cat(3,modelobj.M{simobj.potdep(1:end-1)});
outProj=cat(3,modelobj.outProj{simobj.readouts(2:end)});

updater=mmx('mult',M,outProj);

[alpha,eta,updater]=BWalphaNloop_mex(alpha,eta,updater);


end

