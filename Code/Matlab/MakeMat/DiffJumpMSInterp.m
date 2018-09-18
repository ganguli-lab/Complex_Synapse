function [ Wp,Wm,w ] = DiffJumpMSInterp( n,varargin )
%[Wp,Wm,w]=DIFFJUMP(n) Possibly better than uniform SMS
%   Detailed explanation goes here

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='DiffJumpMSInterp';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addRequired('n',@(x)validateattributes(x,{'numeric'},{'scalar','even'},'DiffJumpMSInterp','n',1))
    p.addOptional('eps',1,@(x)validateattributes(x,{'numeric'},{'scalar','positive'},'DiffJumpMSInterp','eps',2));
end
p.parse(n,varargin{:});
n=p.Results.n;
eps=p.Results.eps;

q=[ones(1,n/2-1) 1-eps ones(1,n/2-1)];

[Wp,Wm,w]=MakeSMS(q);

% Wp(n/2,n/2+1)=0;
Wp(n/2,n)=eps;
% Wm(n/2+1,n/2)=0;
Wm(n/2+1,1)=eps;

Wp=StochastifyC(Wp);
Wm=StochastifyC(Wm);

end

