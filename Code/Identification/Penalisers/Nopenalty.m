function [ newmodelobj,loglike,pstate ] = Nopenalty( modelobj,simobj,varargin )
%[newmodelobj,loglike,pstate]=NOPENALTY(modelobj,simobj) update of
%estiamted HMM with no penalty
%   newmodelobj = updated SynapseIdModel
%   loglike     = log likelihood of readouts under old model (prod over simobj)
%   pstate      = posterior prob of HMM being in each state at each time (cell, one element for each simobj)
%   modelobj = SynapseIdModel
%   simobj   = vector of SynapsePlastSeq

Weighter='RJ';%weighting algorithm to apply before penalty
Normalise=true;
varargin=assignApplicable(varargin);

weighterfn=str2func([Weighter 'weight']);

[newmodelobj,loglike,pstate] = weighterfn(modelobj,simobj,varargin{:});

if Normalise
    newmodelobj=newmodelobj.Normalise;
    assert(newmodelobj.isvalid,'newmodelobj is invalid');
end


end

