function [ newmodelobj,loglike,pstate ] = Mackayweight( modelobj,simobj,varargin )
%[newmodelobj,loglike,pstate]=MACKAYWEIGHT(modelobj,simobj) Mackay update of estimated HMM
%   newmodelobj = updated SynapseIdModel
%   loglike     = log likelihood of readouts under old model (prod over simobj)
%   pstate      = posterior prob of HMM being in each state at each time (cell, one element for each simobj)
%   modelobj = SynapseIdModel
%   simobj   = vector of SynapsePlastSeq

% lloffset=0;%avoid overflow by making this more positive
Algorithm='BW';%algorithm to apply to each chunk of data
HoldBack=-0.25;
Normalise=true;
varargin=assignApplicable(varargin);

updater=str2func([Algorithm 'update']);

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


error(CheckSize(simobj,@isvector));

pstate=cell(1,length(simobj));
newmodelobj=modelobj.Zero;
loglike=0;

for i=1:length(simobj)
    [chunkModel,chunkll,chunkPs]=updater(modelobj,simobj(i),'Normalise',Normalise,varargin{:});
    chunkModel=exp(lloffset+HMMloglike(chunkModel,weightsim))*chunkModel;
    newmodelobj=newmodelobj+chunkModel;
    pstate{i}=chunkPs*diag(1./sum(chunkPs,1));
    loglike=loglike+chunkll;
end
    
newmodelobj=newmodelobj.Normalise;

assert(newmodelobj.isvalid,'newmodelobj is invalid');

end

