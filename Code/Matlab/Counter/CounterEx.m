function [ Wp,Wm,w ] = CounterEx( n,varargin )
%[Wp,Wm,w]=COUNTEREX(n,eps) Counterexample generator
%   Detailed explanation goes here

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='CounterEx';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addRequired('n',@(x)validateattributes(x,{'numeric'},{'scalar','even'},'CounterEx','n',1))
    p.addOptional('eps',1,@(x)validateattributes(x,{'numeric'},{'scalar','positive'},'CounterEx','eps',2));
end
p.parse(n,varargin{:});
r=p.Results;
if any(strcmp('eps',p.UsingDefaults))
    r.eps=1/r.n;
end

q=ones(1,r.n-1);
q(r.n/2+1)=r.eps;
% q(1)=1/n;
[Wp,Wm,w]=MakeSMS(q);


end

