function [ xnew,newfval,exitflag ] = GradientDescend( x,eta,A,b,varargin )
%xnew=GRADIENTDESCEND(x,eta,A,b,fun,extraArgs...) many steps of gradient descent with linear
%constraints
%   x   = params
%   dx  = gradient wrt params
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
    p.FunctionName='GradientDescend';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addRequired('x',@(x)validateattributes(x,{'numeric'},{'column'},'GradientDescend','x',1))
    p.addRequired('eta',@(x)validateattributes(x,{'numeric'},{'scalar'},'GradientDescend','eta',2))
    p.addRequired('A',@(x)validateattributes(x,{'numeric'},{'2d'},'GradientDescend','A',3))
    p.addRequired('b',@(x)validateattributes(x,{'numeric'},{'column'},'GradientDescend','b',4))
    p.addOptional('fun',@OptFunGrad,@(x)validateattributes(x,{'function_handle','char'},{},'GradientDescend','fun',5));
    p.addOptional('extraArgs',{},@(x)validateattributes(x,{'cell'},{},'GradientDescend','extraArgs',5));
    p.addParameter('TolFun',true,@(x) validateattributes(x,{'numeric'},{'scalar'},'GradientDescend','TolFun'));
    p.addParameter('TolX',true,@(x) validateattributes(x,{'numeric'},{'scalar'},'GradientDescend','TolX'));
    p.addParameter('TolCon',true,@(x) validateattributes(x,{'numeric'},{'scalar'},'GradientDescend','TolCon'));
    p.addParameter('MaxIter',true,@(x) validateattributes(x,{'numeric'},{'scalar','integer'},'GradientDescend','MaxIter'));
end
p.parse(x,eta,A,b,varargin{:});
r=p.Results;
error(CheckSize(r.A,@(y)size(y,2)==length(r.x),'samesize(x)'));
error(CheckSize(r.A,@(y)size(y,1)==length(r.b),'samesize(b)'));
if ischar(r.fun) && isrow(r.fun)
    r.fun=str2func(r.fun);
end
x=r.x;

fval=r.fun(x,r.extraArgs{:});
exitflag=0;
for i=1:r.MaxIter
    xnew=GradientStep(x,r.eta,r.A,r.b,r.fun,r.extraArgs,p.Unmatched,'TolCon',r.TolCon);
    newfval=r.fun(xnew,r.extraArgs{:});
    if all(abs(xnew-x))<r.TolX
        exitflag=2;
        break;
    end
    if abs(fval-newfval)<r.TolFun
        exitflag=3;
        break;
    end
    x=xnew;
    fval=newfval;
end



end

