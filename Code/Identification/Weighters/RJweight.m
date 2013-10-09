function [ newmodelobj,loglike,pstate ] = RJweight( modelobj,simobj,varargin )
%[newmodelobj,loglike,pstate]=RJWEIGHT(modelobj,simobj) Rabiner-Juang update of estimated HMM
%   newmodelobj = updated SynapseIdModel
%   loglike     = log likelihood of readouts under old model (prod over simobj)
%   pstate      = posterior prob of HMM being in each state at each time (cell, one element for each simobj)
%   modelobj = SynapseIdModel
%   simobj   = vector of SynapsePlastSeq

% lloffset=0;%avoid overflow by making this more positive
Algorithm='BW';%algorithm to apply to each chunk of data
Normalise=false;
varargin=assignApplicable(varargin);

updater=str2func([Algorithm 'update']);

error(CheckSize(simobj,@isvector));

lloffset=-HMMloglike(modelobj,simobj)/length(simobj);

pstate=cell(1,length(simobj));
newmodelobj=modelobj.Zero;
loglike=0;

for i=1:length(simobj)
    [chunkModel,chunkll,chunkPs]=updater(modelobj,simobj(i),'Normalise',Normalise,varargin{:});
    if ~strcmpi(Algorithm,'BW')
        chunkModel=exp(-lloffset-chunkll)*chunkModel;
    end
    newmodelobj=newmodelobj+chunkModel;
    pstate{i}=chunkPs*diag(1./sum(chunkPs,1));
    loglike=loglike+chunkll;
end
    
newmodelobj=newmodelobj.Normalise;

assert(newmodelobj.isvalid,'newmodelobj is invalid');

end

