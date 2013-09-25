function [ newmodelobj,pstate,loglike ] = Uniweight( modelobj,simobj,varargin )
%[newmodelobj,pstate,loglike]=UNIWEIGHT(modelobj,simobj) update of
%estiamted HMM with equal weights
%   newmodelobj = updated SynapseIdModel
%   pstate      = posterior prob of HMM being in each state at each time (cell, one element for each simobj)
%   loglike     = log likelihood of readouts under old model (prod over simobj)
%   modelobj = SynapseIdModel
%   simobj   = vector of SynapsePlastSeq

Algorithm='BW';%algorithm to apply to each chunk of data
varargin=assignApplicable(varargin);

updater=str2func([Algorithm 'update']);

error(CheckSize(simobj,@isvector));

pstate=cell(1,length(simobj));
newmodelobj=modelobj.Zero;
loglike=0;

for i=1:length(simobj)
    [chunkModel,chunkPs,chunkll]=updater(modelobj,simobj(i),varargin{:});
    newmodelobj=newmodelobj+chunkModel;
    pstate{i}=chunkPs*diag(1./sum(chunkPs,1));
    loglike=loglike+chunkll;
end
    
newmodelobj=newmodelobj.Normalise;

assert(newmodelobj.isvalid,'newmodelobj is invalid');

end

