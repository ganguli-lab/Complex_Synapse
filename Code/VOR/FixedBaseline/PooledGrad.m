function [ dWp,dWm ] = PooledGrad( q,n )
%[dWp,dWm]=POOLEDGRAD(q,n)  gradient of matrix elements wrt parameter
%   q  = ratio of nearby transition rates
%   n  = number of states (even)
%   dWP = gradient of potentiation transition rates
%   dWM = gradient of depression transition rates


[dWp,dWm] = SerialBuilder(n,1);

end

