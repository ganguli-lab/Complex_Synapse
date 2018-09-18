function [ pinf ] = EqProb( W,varargin )
%EQPROB(W) equlibrium distribution of Markov chain (cts time)
%   W = transition rates

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='EqProb';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addRequired('W',@(x)validateattributes(x,{'numeric'},{'2d','square'},'EqProb','W',1));
    p.addParameter('RCondThresh',1e-5,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'},'EqProb','RCondThresh'));
end
p.parse(W,varargin{:});
r=p.Results;


Zinv=ones(length(r.W)) - r.W;

if rcond(Zinv)>r.RCondThresh
    pinf = ones(1,size(r.W,1))/(Zinv);
else
    [v,qb]=eig(-r.W');
    qb=diag(qb);
    [~,ix]=min(qb);
    pinf=v(:,ix).';
    pinf=pinf/sum(pinf);
end

end

