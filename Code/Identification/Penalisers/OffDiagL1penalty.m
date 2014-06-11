function [ newmodelobj,loglike,pstate ] = OffDiagL1penalty( modelobj,simobj,varargin )
%[newmodelobj,loglike,pstate]=OFFDIAGL1PENALTY(modelobj,simobj) update of
%estiamted HMM with off diagonal L^1 penalty
%   newmodelobj = updated SynapseIdModel
%   loglike     = log likelihood of readouts under old model (prod over simobj)
%   pstate      = posterior prob of HMM being in each state at each time (cell, one element for each simobj)
%   modelobj = SynapseIdModel
%   simobj   = vector of SynapsePlastSeq

Weighter='RJ';%weighting algorithm to apply before penalty
Normalise=true;
Penalty=1;
varargin=assignApplicable(varargin);

weighterfn=str2func([Weighter 'weight']);

[newmodelobj,loglike,pstate] = weighterfn(modelobj,simobj,'Normalise',false,varargin{:});

M=newmodelobj.M;
L1=0;
for i=1:newmodelobj.NumPlast
    lambda=CalcLagrange(M{i});
    M{i}=PenNorm(M{i},lambda);
    L1=L1+sum(sum(M{i}))-trace(M{i});
end
newmodelobj=newmodelobj.setM(M);

loglike = loglike - Penalty * L1;

if Normalise
    newmodelobj=newmodelobj.Normalise;
    assert(newmodelobj.isvalid,'newmodelobj is invalid');
end

    function lambda=CalcLagrange(M)
        rowsum = sum(M,2) - Penalty;
        lambda = ( rowsum + sqrt( rowsum.^2 + 4*Penalty*diag(M) ) )/2;
    end

    function Mnew=PenNorm(M,lambda)
        ev=ones(wrev(size(lambda)));
        Mnew = M ./ ( (lambda+Penalty) * ev );
        Mnew = Mnew - diag(diag(Mnew)) + diag( diag(M)./lambda );
    end

end

