function [ pinf ] = EqProb( obj,varargin )
%pinf=SynapseMemoryModel.EQPROB equlibrium distribution of Markov chain (cts time)


RCondThresh=1e-3;
varargin=assignApplicable(varargin);

[Zinv,piv]=obj.GetZinv;

if rcond(Zinv)>RCondThresh
    pinf = piv/Zinv;
else
    [v,qb]=eig(-obj.GetWf');
    qb=diag(qb);
    [~,ix]=min(qb);
    pinf=v(:,ix).';
    pinf=pinf/sum(pinf);
end

end

