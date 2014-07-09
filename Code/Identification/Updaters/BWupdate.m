function [ newmodelobj,loglike,pstate ] = BWupdate( modelobj,simobj,varargin )
%[newmodelobj,loglike,pstate]=BWUPDATE(modelobj,simobj) Baum-Welch update of estimated HMM
%   newmodelobj = updated SynapseIdModel
%   loglike     = log likelihood of readout given current model
%   pstate      = posterior prob of HMM being in each state at each time
%   modelobj = SynapseIdModel
%   simobj   = SynapsePlastSeq



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


[alpha,eta,updater]=BWalphaN(modelobj,simobj);
beta=BWbetaN(updater,modelobj,simobj);

pstate=alpha'.*beta;

alpha=reshape(alpha(1:end-1,:)', [size(alpha,2) 1 size(alpha,1)-1]);
beta=reshape(beta(:,2:end), [1 size(beta,1) size(beta,2)-1]);

Mupdate = mmx('mult',alpha,beta) .* updater;

M_new=cell(1,modelobj.NumPlast);
for i=1:length(M_new)
    M_new{i} =  sum( Mupdate(:,:,simobj.potdep(1:end-1)==i) , 3 )  ;
end
    


newmodelobj=modelobj.setM(M_new);
newmodelobj=newmodelobj.setInitial(pstate(:,1)');

if p.Results.Normalise
    pstate=pstate*diag(1./sum(pstate,1));
    newmodelobj=newmodelobj.Normalise;
    assert(newmodelobj.isvalid,'newmodelobj is invalid');
end


loglike = -sum(log(eta));


end

