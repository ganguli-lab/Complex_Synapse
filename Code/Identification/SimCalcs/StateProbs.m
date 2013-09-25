function [ p ] = StateProbs( modelobj,simobj )
%p=STATEPROBS(modelobj,simobj) posterior prob of HMM being in each state at each time
%   modelobj = SynapseIdModel
%   simobj   = SynapsePlastSeq



p=zeros(length(modelobj.Initial),length(simobj.readouts));
% alpha=BWalpha(length(readouts),readouts,Initial,outProj,M,potdep);
% beta=BWbeta(1,readouts,outProj,M,potdep);
[alpha,eta]=BWalphaN(length(simobj.readouts),modelobj,simobj);
beta=BWbetaN(eta,1,modelobj,simobj);

for t=1:length(simobj.readouts)
    p(:,t)=alpha(t,:)'.*beta(:,t);
    p(:,t)=p(:,t)/sum(p(:,t));
end

end

