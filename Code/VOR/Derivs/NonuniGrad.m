function [ dWp, dWm ] = NonuniGrad( n, x )
%[dWp,dWm]=NONUNIGRAD(n,x) gradient of matrix elements wrt parameter
%   n  = number of states (even)
%   x  = ratio of nearby transition rates
%   dWP = gradient of potentiation transition rates
%   dWM = gradient of depression transition rates

expnt = abs((1:n-1)-n/2);
% qp=x.^expnt;
dqp = expnt .* x.^(expnt - 1);

[dWp,dWm]=MakeSMS(dqp);


end

