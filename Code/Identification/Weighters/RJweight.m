function [ newmodelobj,pstate,loglike ] = RJweight( modelobj,simobj,varargin )
%[newmodelobj,pstate,loglike]=RJUPDATE(modelobj,simobj) Rabiner-Juang update of estiamted HMM
%   newmodelobj = updated SynapseIdModel
%   pstate      = posterior prob of HMM being in each state at each time (cell, one element for each simobj)
%   loglike     = log likelihood of readouts under old model (prod over simobj)
%   modelobj = SynapseIdModel
%   simobj   = vector of SynapsePlastSeq

lloffset=0;%avoid overflow by making this more positive
algorithm='BW';%algorithm to apply to each chunk of data
varargin=assignApplicable(varargin);

updater=str2func([algorithm 'update']);

error(CheckSize(simobj,@isvector));

pstate=cell(1,length(simobj));
newmodelobj=modelobj.Zero;
loglike=0;

for i=1:length(simobj)
    [chunkModel,chunkPs,chunkll]=updater(modelobj,simobj(i),'Normalise',false,varargin{:});
    Weight=exp(-lloffset-chunkll);
    newmodelobj=newmodelobj+Weight*chunkModel;
    pstate{i}=chunkPs*diag(1./sum(chunkPs,1));
    loglike=loglike+chunkll;
end
    
newmodelobj=newmodelobj.Normalise;

assert(newmodelobj.isvalid,'newmodelobj is invalid');

end

