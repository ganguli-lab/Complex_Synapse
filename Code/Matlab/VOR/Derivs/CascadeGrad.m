function [ dWp, dWm ] = CascadeGrad( n, x )
%[dWp,dWm]=CASCGRAD(n,x) gradient of matrix elements wrt parameter
%   n  = number of states (even)
%   x  = ratio of nearby transition rates
%   dWP = gradient of potentiation transition rates
%   dWM = gradient of depression transition rates

%qp=(xp.^abs((1:n-1)-n/2))/(1-xp);
expnt = abs((1:n-1)-n/2);
dqp = expnt .* (x.^(expnt - 1)) / (1-x) + (x.^expnt) / (1 - x)^2;
[dWp,dWm]=CascadeMSinterp(dqp,dqp,0.5,1);


end

