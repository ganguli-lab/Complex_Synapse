function [ dWp,dWm ] = PooledGrad( n,q )
%[dWp,dWm]=POOLEDGRAD(n,q)  gradient of matrix elements wrt parameter
%   n  = number of states
%   q  = range of transition rates
%   dWP = gradient of potentiation transition rates
%   dWM = gradient of depression transition rates

qq = (n-1:-1:1)/(n-1);

if isscalar(q)
    [dWp,dWm] = MakeSMS(qq);
else
    [dpmax,dmmax] = MakeSMS(qq.^2);
    [dpmin,dmmin] = MakeSMS(qq .* wrev(qq));
    dWp = {dpmax, dpmin};
    dWm = {dmmax, dmmin};
end

end

