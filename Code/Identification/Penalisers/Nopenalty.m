function [ newmodelobj,loglike,pstate ] = Nopenalty( modelobj,simobj,varargin )
%[newmodelobj,loglike,pstate]=NOPENALTY(modelobj,simobj) update of
%estiamted HMM with no penalty
%   newmodelobj = updated SynapseIdModel
%   loglike     = log likelihood of readouts under old model (prod over simobj)
%   pstate      = posterior prob of HMM being in each state at each time (cell, one element for each simobj)
%   modelobj = SynapseIdModel
%   simobj   = vector of SynapsePlastSeq

% Weighter='RJ';%weighting algorithm to apply before penalty
% Normalise=true;
% varargin=assignApplicable(varargin);
persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='Nopenalty';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addParameter('Weighter','RJ');
    p.addParameter('Normalise',true);
    p.addParameter('Penalty',1);
%     p.addParameter('Normalise',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
end
p.parse(varargin{:});

weighterfn=str2func([p.Results.Weighter 'weight']);

[newmodelobj,loglike,pstate] = weighterfn(modelobj,simobj,'Normalise',false,p.Unmatched);

if p.Results.Normalise
    newmodelobj=newmodelobj.Normalise;
%     assert(newmodelobj.isvalid,'newmodelobj is invalid');
end


end

