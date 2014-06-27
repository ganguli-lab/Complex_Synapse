function [ newmodelobj,loglike,pstate ] = RJweight( modelobj,simobj,varargin )
%[newmodelobj,loglike,pstate]=RJWEIGHT(modelobj,simobj) Rabiner-Juang update of estimated HMM
%   newmodelobj = updated SynapseIdModel
%   loglike     = log likelihood of readouts under old model (prod over simobj)
%   pstate      = posterior prob of HMM being in each state at each time (cell, one element for each simobj)
%   modelobj = SynapseIdModel
%   simobj   = vector of SynapsePlastSeq

% lloffset=0;%avoid overflow by making this more positive
% Algorithm='BW';%algorithm to apply to each chunk of data
% Normalise=true;
% varargin=assignApplicable(varargin);
persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='RJweight';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addParameter('Algorithm','BW');
    p.addParameter('Normalise',true);
%     p.addParameter('Normalise',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
end
p.parse(varargin{:});

updater=str2func([p.Results.Algorithm 'update']);

% error(CheckSize(simobj,@isvector));

notBW=~strcmpi(p.Results.Algorithm(end-1:end),'BW');

if notBW
    lloffset=-HMMloglike(modelobj,simobj)/length(simobj);
end

pstate=cell(1,length(simobj));
newmodelobj=modelobj.Zero;
loglike=0;

for i=1:length(simobj)
    if notBW
        [chunkModel,chunkll,chunkPs]=updater(modelobj,simobj(i),'Normalise',true,p.Unmatched);
        chunkModel=exp(-lloffset-chunkll)*chunkModel;
    else
        [chunkModel,chunkll,chunkPs]=updater(modelobj,simobj(i),'Normalise',false,p.Unmatched);
    end
    newmodelobj=newmodelobj+chunkModel;
    pstate{i}=chunkPs*diag(1./sum(chunkPs,1));
    loglike=loglike+chunkll;
end
    
if p.Results.Normalise
    newmodelobj=newmodelobj.Normalise;
    assert(newmodelobj.isvalid,'newmodelobj is invalid');
end

end

