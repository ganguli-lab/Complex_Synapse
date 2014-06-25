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
% siz=modelobj.NumStates*[1 1];
% for i=2:simobj.NumT
%     alpha(i,:) = alpha(i-1,:) * reshape(updater(:,:,i-1),siz);
%     eta(i) = 1 / sum(alpha(i,:));
%     alpha(i,:) = alpha(i,:) * eta(i);
% end
[ alpha,eta,updater ] = BWalphaNloop_mex( alpha,eta,updater );
% updater = updater .* reshape( ones((modelobj.NumStates )^2,1) * eta(2:end)', [modelobj.NumStates modelobj.NumStates simobj.NumT-1] );


end

