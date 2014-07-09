function [ newmodelobj,loglike,pstate ] = Mackayweight( modelobj,simobj,varargin )
%[newmodelobj,loglike,pstate]=MACKAYWEIGHT(modelobj,simobj) Mackay update of estimated HMM
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
    p.addParameter('HoldBack',-0.25);
    p.addParameter('Normalise',true);
%     p.addParameter('Normalise',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
end
p.parse(varargin{:});

updater=str2func([p.Results.Algorithm 'update']);

HoldBack=p.Results.HoldBack;
if isscalar(HoldBack) && HoldBack <0
    HoldBack=-HoldBack;
    if ~isint(HoldBack)
        HoldBack=ceil(HoldBack*length(simobj));
    end%isint
    HoldBack=randperm(length(simobj),HoldBack);
end

weightsim=simobj(HoldBack);
simobj(HoldBack)=[];

lloffset=-HMMloglike(modelobj,weightsim)/length(weightsim);


% error(CheckSize(simobj,@isvector));

pstate=cell(1,length(simobj));
newmodelobj=modelobj.Zero;
loglike=0;

for i=1:length(simobj)
    [chunkModel,chunkll,chunkPs]=updater(modelobj,simobj(i),'Normalise',true,p.Unmatched);
    chunkModel=exp(lloffset+HMMloglike(chunkModel,weightsim))*chunkModel;
    newmodelobj=newmodelobj+chunkModel;
    pstate{i}=chunkPs*diag(1./sum(chunkPs,1));
    loglike=loglike+chunkll;
end
    
if p.Results.Normalise
    newmodelobj=newmodelobj.Normalise;
    assert(newmodelobj.isvalid,'newmodelobj is invalid');
end

end

