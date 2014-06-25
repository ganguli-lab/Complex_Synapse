function [ newmodelobj,loglike,pstate ] = InitBWupdate( modelobj,simobj,varargin )
%[newmodelobj,loglike,pstate]=INITBWUPDATE(modelobj,simobj) Baum-Welch update of estimated HMM
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
    p.FunctionName='InitBWupdate';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addParameter('Normalise',true);
%     p.addParameter('Normalise',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
end
p.parse(varargin{:});


% pstate=zeros(length(modelobj.Initial),length(simobj.readouts));
[alpha,eta,updater]=BWalphaN(modelobj,simobj);
beta=BWbetaN(updater,modelobj,simobj);


% for t=1:length(simobj.readouts)
%     pstate(:,t)=alpha(t,:)'.*beta(:,t);
% end
pstate=alpha'.*beta;

newmodelobj=modelobj.setInitial(pstate(:,1)');

if p.Results.Normalise
    pstate=pstate*diag(1./sum(pstate,1));
    newmodelobj=newmodelobj.Normalise;
    assert(newmodelobj.isvalid,'newmodelobj is invalid');
end


loglike = -sum(log(eta));


end

