function [ beta ] = BWbetaN( updater,modelobj,simobj )
%beta=BWBETA(upsater,modelobj,simobj) normalised backward variables for Baum-Welch algorithm
%   beta     = normalised backward variables
%   updater  = updater matrices (M*outProj*eta)
%   modelobj = SynapseIdModel
%   simobj   = SynapsePlastSeq

%   t        = time-step after which we want forward variables

% error(CheckSize(modelobj,@isscalar));
% error(CheckSize(simobj,@isscalar));

% error(CheckSize(t,@isscalar));
% error(CheckValue(t,@isint));

% M=cat(3,modelobj.M{simobj.potdep});
% outProj=cat(3,modelobj.outProj{simobj.readouts});

beta=ones(modelobj.NumStates,simobj.NumT);

%     siz=modelobj.NumStates*[1 1];
%     for i=simobj.NumT:-1:2
%         beta(:,i-1) = reshape(updater(:,:,i-1),siz) * beta(:,i);
%     end
[ beta ] = BWbetaNloop_mex( beta,updater );

end

