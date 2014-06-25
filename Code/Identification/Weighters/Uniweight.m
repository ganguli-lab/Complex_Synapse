function [ newmodelobj,loglike,pstate ] = Uniweight( modelobj,simobj,varargin )
%[newmodelobj,loglike,pstate]=UNIWEIGHT(modelobj,simobj) update of
%estiamted HMM with equal weights
%   newmodelobj = updated SynapseIdModel
%   loglike     = log likelihood of readouts under old model (prod over simobj)
%   pstate      = posterior prob of HMM being in each state at each time (cell, one element for each simobj)
%   modelobj = SynapseIdModel
%   simobj   = vector of SynapsePlastSeq

% Algorithm='BW';%algorithm to apply to each chunk of data
% Normalise=true;
% varargin=assignApplicable(varargin);
persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='Uniweight';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addParameter('Algorithm','BW');
    p.addParameter('Normalise',true);
%     p.addParameter('Normalise',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
end
p.parse(varargin{:});

updater=str2func([p.Results.Algorithm 'update']);

% error(CheckSize(simobj,@isvector));

pstate=cell(1,length(simobj));
newmodelobj=modelobj.Zero;
loglike=0;

for i=1:length(simobj)
    [chunkModel,chunkll,chunkPs]=updater(modelobj,simobj(i),p.Unmatched);
    newmodelobj=newmodelobj+chunkModel;
    pstate{i}=chunkPs*diag(1./sum(chunkPs,1));
    loglike=loglike+chunkll;
end
    
if p.Results.Normalise
    newmodelobj=newmodelobj.Normalise;
    assert(newmodelobj.isvalid,'newmodelobj is invalid');
end

end

