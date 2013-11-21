function [ xnew ] = GradientStep( x,eta,A,b,fun,varargin )
%xnew=GRADIENTSTEP(x,eta,A,b,fun) one step of gradient descent with linear
%constraints
%   x   = params
%   fun  = gradient wrt params
%   eta = learning rate
% dx/dt = -eta * (dx + (mu*A)')
%   A   = gradient of constraints
%   b   = intercept of constraints
% A * x <= b
% L = fun(x) + mu*(A*x-b)
% dx = dfun(x)/dx

error(CheckSize(x,@iscolumn));
error(CheckSize(eta,@isscalar));
error(CheckSize(A,@ismat));
error(CheckSize(b,@iscolumn));
error(CheckSize(A,@(y)size(y,2)==length(x),'samesize(x)'));
error(CheckSize(A,@(y)size(y,1)==length(b),'samesize(b)'));

existsAndDefault('fun',@OptFunGrad);
if ischar(fun) && isrow(fun)
    fun=str2func(fun);
end
error(CheckType(fun,'function_handle'));

TolCon = 1e-10;
varargin=assignApplicable(varargin);




[~,dx]=fun(x,varargin{:});
%dx=-dx;
mu=CalcKTmult(x,dx,eta,A,b);

xnew = x - eta *(dx + (mu*A).');


end

