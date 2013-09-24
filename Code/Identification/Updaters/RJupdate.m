function [ newmodelobj,pstate,loglike ] = RJupdate( modelobj,simobj,chunks,varargin )
%[M_new,initial_new,pstate,loglike]=RJUPDATE(chunks,readouts,initial,outProj,M,potdep)
%Rabiner-Juang update of estiamted HMM
%   newmodelobj = updated SynapseIdModel
%   pstate      = posterior prob of HMM being in each state at each time
%   loglike     = log likelihood of readouts under old model (prod over chunks)
%   chunks   = 2-by-K matrix of starts and ends of each chunk.
%   modelobj = SynapseIdModel
%   simobj   = SynapsePlastSeq

lloffset=0;%avoid overflow by making this more positive
varargin=assignApplicable(varargin);


pstate=zeros(length(modelobj.initial),length(simobj.readouts));
newmodelobj=modelobj.Zero;
loglike=0;

for i=1:size(chunks,2)
    [chunkModel,chunkPs,chunkll]=BWupdate(modelobj,simobj.GetRange(chunks(:,i)),'Normalise',false);
    Weight=exp(-lloffset-chunkll);
    newmodelobj=newmodelobj+Weight*chunkModel;
    pstate(:,chunks(1,i):chunks(2,i))=chunkPs*diag(1./sum(chunkPs,1));
    loglike=loglike+chunkll;
end
    
newmodelobj=newmodelobj.Normalise;

assert(newmodelobj.isvalid,'newmodelobj is invalid');

end

