function [ newmodelobj,loglike,pstate ] = Sequentialweight( modelobj,simobj,varargin )
%[newmodelobj,loglike,pstate]=SEQUENTIAlWEIGHT(modelobj,simobj) Sequential update of estimated HMM
%   newmodelobj = updated SynapseIdModel
%   loglike     = log likelihood of readouts under old model (prod over simobj)
%   pstate      = posterior prob of HMM being in each state at each time (cell, one element for each simobj)
%   modelobj = SynapseIdModel
%   simobj   = vector of SynapsePlastSeq

% lloffset=0;%avoid overflow by making this more positive
persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='Mackayweight';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addParameter('Algorithm','BW');
    p.addParameter('RecurrentFrac',0.1);
    p.addParameter('Normalise',true);
%     p.addParameter('Normalise',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
end
p.parse(varargin{:});

updater=str2func([p.Results.Algorithm 'update']);



% error(CheckSize(simobj,@isvector));

pstate=cell(1,length(simobj));
newmodelobj=modelobj;
loglike=0;

for i=1:length(simobj)
    [chunkModel,chunkll,chunkPs]=updater(newmodelobj,simobj(i),'Normalise',false,p.Unmatched);
    newmodelobj=newmodelobj+chunkModel +...
        max(cellfun(@(x)max(max(x)),chunkModel.M))*p.Results.RecurrentFrac*SynapseIdModel.Rand(newmodelobj.w);
    newmodelobj=newmodelobj.Normalise;
    pstate{i}=chunkPs*diag(1./sum(chunkPs,1));
    loglike=loglike+chunkll;
end
    
if p.Results.Normalise
    newmodelobj=newmodelobj.Normalise;
    assert(newmodelobj.isvalid,'newmodelobj is invalid');
end

end

