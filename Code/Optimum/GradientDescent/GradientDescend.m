function [ xnew,newfval,exitflag ] = GradientDescend( x,eta,A,b,fun,varargin )
%xnew=GRADIENTDESCEND(x,eta,A,b,fun,...) many steps of gradient descent with linear
%constraints
%   x   = params
%   dx  = gradient wrt params
%   eta = learning rate
% dx/dt = -eta * (dx + (mu*A)')
%   A   = gradient of constraints
%   b   = intercept of constraints
% A * x <= b
% L = fun(x) + mu*(A*x-b)
% dx = dfun(x)/dx

error(CheckSize(x,@iscolumn));
error(CheckSize(eta,@isscalar));

existsAndDefault('fun',@OptFunGrad);
if ischar(fun) && isrow(fun)
    fun=str2func(fun);
end
error(CheckType(fun,'function_handle'));

TolFun=1e-6;
TolX=1e-10;
TolCon=1e-6;
MaxIter=1000;
varargin=assignApplicable(varargin);

fval=fun(x,varargin{:});
exitflag=0;
for i=1:MaxIter
    xnew=GradientStep(x,eta,A,b,fun,'TolCon',TolCon,varargin{:});
    newfval=fun(xnew,varargin{:});
    if all(abs(xnew-x))<TolX
        exitflag=2;
        break;
    end
    if abs(fval-newfval)<TolFun
        exitflag=3;
        break;
    end
    x=xnew;
    fval=newfval;
end



end

