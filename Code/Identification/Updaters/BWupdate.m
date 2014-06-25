function [ newmodelobj,loglike,pstate ] = BWupdate( modelobj,simobj,varargin )
%[newmodelobj,loglike,pstate]=BWUPDATE(modelobj,simobj) Baum-Welch update of estimated HMM
%   newmodelobj = updated SynapseIdModel
%   loglike     = log likelihood of readout given current model
%   pstate      = posterior prob of HMM being in each state at each time
%   modelobj = SynapseIdModel
%   simobj   = SynapsePlastSeq


% Normalise=true;
% varargin=assignApplicable(varargin);

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='BWupdate';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addParameter('Normalise',true);
%     p.addParameter('Normalise',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
end
p.parse(varargin{:});


siz=modelobj.NumStates*[1 1];
[alpha,eta,updater]=BWalphaN(modelobj,simobj);
beta=BWbetaN(updater,modelobj,simobj);
% [alpha,eta]=BWalphaN(modelobj,simobj);
% beta=BWbetaN(eta,modelobj,simobj);
% M_new={zeros(length(modelobj.M{1}))};
% for i=2:numel(modelobj.M)
%     M_new{i}=M_new{1};
% end


% for t=1:simobj.NumT-1
%     M_new{simobj.potdep(t)}=M_new{simobj.potdep(t)} + (beta(:,t+1)*alpha(t,:))' .*...
%         reshape(updater(:,:,t),siz);
% %     M_new{simobj.potdep(t)}=M_new{simobj.potdep(t)} + (beta(:,t+1)*alpha(t,:))' .*...
% %         (modelobj.M{simobj.potdep(t)}*modelobj.outProj{simobj.readouts(t+1)}) * eta(t+1);
% end
pstate=alpha'.*beta;

alpha=reshape(alpha(1:end-1,:)', [modelobj.NumStates 1 simobj.NumT-1]);
beta=reshape(beta(:,2:end), [1 modelobj.NumStates simobj.NumT-1]);

Mupdate = mmx('mult',alpha,beta) .* updater;
M_new=cell(1,modelobj.NumPlast);
for i=1:modelobj.NumPlast
    M_new{i} = reshape( sum( Mupdate(:,:,simobj.potdep(1:end-1)==i) , 3 )  , siz );
end
    


% newmodelobj=SynapseIdModel(modelobj,'M',M_new,'Initial',pstate(:,1)');
newmodelobj=modelobj.setM(M_new);
newmodelobj=newmodelobj.setInitial(pstate(:,1)');

if p.Results.Normalise
    pstate=pstate*diag(1./sum(pstate,1));
    newmodelobj=newmodelobj.Normalise;
    assert(newmodelobj.isvalid,'newmodelobj is invalid');
end


loglike = -sum(log(eta));


end

