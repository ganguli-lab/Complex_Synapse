function [ Wp,Wm,w ] = DiffJump( n,varargin )
%[Wp,Wm,w]=DIFFJUMP(n) Possibly better than uniform SMS
%   Detailed explanation goes here

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='DiffJump';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addRequired('n',@(x)validateattributes(x,{'numeric'},{'scalar','even'},'DiffJump','W',1))
    p.addOptional('eps',1,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1},'DiffJump','eps',2));
    p.addOptional('q1',1,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1},'DiffJump','q1',3));
    p.addOptional('q2',1,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1},'DiffJump','q2',4));
    p.addOptional('q3',1,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1},'DiffJump','q3',5));
end
p.parse(n,varargin{:});
r=p.Results;
n=r.n;

q=[r.q1*ones(1,n/2-1) 0 r.q2*ones(1,n/2-1)];
q(1)=r.eps;

[Wp,Wm,w]=MakeSMS(q);

% Wp(n/2,n/2+1)=0;
Wp(n/2,n)=r.q3;
% Wm(n/2+1,n/2)=0;
Wm(n/2+1,1)=r.q3;

Wp=StochastifyC(Wp);
Wm=StochastifyC(Wm);

end

