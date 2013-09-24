function [ newmodelobj,pstate,loglike ] = BWupdate( modelobj,simobj,varargin )
%[M_new,initial_new,pstate,loglike]=BWUPDATE(readouts,initial,outProj,M,potdep)
%Baum-Welch update of estiamted HMM
%   newmodelobj = updated SynapseIdModel
%   pstate      = posterior prob of HMM being in each state at each time
%   loglike     = log likelihood of readout given current model
%   modelobj = SynapseIdModel
%   simobj   = SynapsePlastSeq


Normalise=true;
varargin=assignApplicable(varargin);


pstate=zeros(length(modelobj.initial),length(simobj.readouts));
[alpha,eta]=BWalphaN(length(simobj.readouts),modelobj,simobj);
beta=BWbetaN(eta,1,modelobj,simobj);
M_new={zeros(length(modelobj.M{1}))};
for i=2:numel(modelobj.M)
    M_new{i}=M_new{1};
end


for t=1:length(readouts)-1
    pstate(:,t)=alpha(t,:)'.*beta(:,t);
    M_new{simobj.potdep(t)}=M_new{simobj.potdep(t)} + (beta(:,t+1)*alpha(t,:))' .*...
        (simobj.M{simobj.potdep(t)}*simobj.outProj{simobj.readouts(t+1)}) * eta(t+1);
end
    pstate(:,end)=alpha(end,:)'.*beta(:,end);

newmodelobj=SynapseIdModel(modelobj,'M',M_new,'initial',pstate(:,1)');

if Normalise
    pstate=pstate*diag(1./sum(pstate,1));
    newmodelobj=newmodelobj.Normalise;
    assert(newmodelobj.isvalid,'newmodelobj is invalid');
end


loglike = -sum(log(eta));


end
