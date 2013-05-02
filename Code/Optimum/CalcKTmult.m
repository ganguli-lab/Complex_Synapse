function [ mu ] = CalcKTmult( x,dx,eta,A,b )
%mu=CALCKTMULT(x,dx,eta,A,b) Calculate Kuhn-Tucker multipliers for gradient
%descent with linear constraints
%   x   = params
%   dx  = gradient wrt params
%   eta = learning rate
% dx/dt = -eta * (dx + (mu*A)')
%   A   = gradient of constraints
%   b   = intercept of constraints
% A * x <= b
% L = f(x) + mu*(A*x-b)

error(CheckSize(x,@iscolumn));
error(CheckSize(dx,@iscolumn));
error(CheckSize(dx,@(y)samesize(x,y),'samesize(x)'));
error(CheckSize(eta,@isscalar));
error(CheckSize(A,@ismat));
error(CheckSize(b,@iscolumn));
error(CheckSize(A,@(y)size(y,2)==length(x),'samesize(x)'));
error(CheckSize(A,@(y)size(y,1)==length(b),'samesize(b)'));

mu = zeros(1,size(A,1));

ineq = A*(x-eta*dx) - b;
nonslack = (ineq > 0);

An=A(nonslack,:);

mu(nonslack) = (An*An.')\ineq(nonslack)/eta;




end

