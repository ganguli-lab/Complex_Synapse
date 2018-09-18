function [ pinf ] = EqProb( obj,varargin )
%pinf=SynapseMemoryModel.EQPROB equlibrium distribution of Markov chain (cts time)

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='SynapseMemoryModel.EqProb';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addParameter('RCondThresh',1e-5,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'},'SynapseMemoryModel.EqProb','RCondThresh'));
end
p.parse(varargin{:});
r=p.Results;


[Zinv,piv]=obj.GetZinv;

if rcond(Zinv)>r.RCondThresh
    pinf = piv/Zinv;
else
    [v,qb]=eig(-obj.GetWf');
    qb=diag(qb);
    [~,ix]=min(qb);
    pinf=v(:,ix).';
    pinf=pinf/sum(pinf);
end

end

