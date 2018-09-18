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


beta=ones(modelobj.NumStates,simobj.NumT);

beta=BWbetaNloop_mex(beta,updater);

end

