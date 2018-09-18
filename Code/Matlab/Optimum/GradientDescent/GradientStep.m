function [ xnew ] = GradientStep( x,eta,A,b,varargin )
%xnew=GRADIENTSTEP(x,eta,A,b,fun,extraArgs) one step of gradient descent with linear
%constraints
%   x   = params
%   fun  = gradient wrt params
%   eta = learning rate
% dx/dt = -eta * (dx + (mu*A)')
%   A   = gradient of constraints
%   b   = intercept of constraints
%   extraArgs = cell of extra arguments passed to fun
% A * x <= b
% L = fun(x) + mu*(A*x-b)
% dx = dfun(x)/dx

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='GradientStep';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addRequired('x',@(x)validateattributes(x,{'numeric'},{'column'},'GradientStep','x',1))
    p.addRequired('eta',@(x)validateattributes(x,{'numeric'},{'scalar'},'GradientStep','eta',2))
    p.addRequired('A',@(x)validateattributes(x,{'numeric'},{'2d'},'GradientStep','A',3))
    p.addRequired('b',@(x)validateattributes(x,{'numeric'},{'column'},'GradientStep','b',4))
    p.addOptional('fun',@OptFunGrad,@(x)validateattributes(x,{'function_handle','char'},{},'GradientStep','fun',5));
    p.addOptional('extraArgs',{},@(x)validateattributes(x,{'cell'},{},'GradientStep','extraArgs',5));
    p.addParameter('TolCon',1e-10,@(x) validateattributes(x,{'numeric'},{'scalar'},'GradientStep','TolCon'));
end
p.parse(x,eta,A,b,varargin{:});
r=p.Results;
error(CheckSize(r.A,@(y)size(y,2)==length(r.x),'samesize(x)'));
error(CheckSize(r.A,@(y)size(y,1)==length(r.b),'samesize(b)'));
if ischar(r.fun) && isrow(r.fun)
    r.fun=str2func(r.fun);
end
x=r.x;




[~,dx]=r.fun(x,r.extraArgs{:});
%dx=-dx;
mu=CalcKTmult(x,dx,r.eta,r.A,r.b);

xnew = x - r.eta *(dx + (mu*r.A).');


end

