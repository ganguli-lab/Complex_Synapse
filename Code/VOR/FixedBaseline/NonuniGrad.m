function [ dWp, dWm ] = NonuniGrad( x,n )
%[dWp,dWm]=NONUNIGRAD(x,n) gradient of matrix elements wrt parameter
%   x  = ratio of nearby transition rates
%   n  = number of states (even)
%   dWP = gradient of potentiation transition rates
%   dWM = gradient of depression transition rates

expnt = abs((1:n-1)-n/2);
% qp=x.^expnt;
dqp = expnt .* x.^(expnt - 1);

[dWp,dWm]=MakeSMS(dqp);


end

