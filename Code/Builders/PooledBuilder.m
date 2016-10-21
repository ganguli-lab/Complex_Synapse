function [ Wp,Wm,w ] = PooledBuilder( n,q )
%[Wp,Wm,w]=POOLEDBUILDER(n,q) Build pooled resource model
%   q  = range of transition rates [qmax qmin]
%   n  = number of states (even)
%   Wp = potentiation transition rates
%   Wm = depression transition rates
%   w  = Weights of states (+/-1)

    w=(-1:2/(n-1):1)';
    if isscalar(q) || diff(q)==0
        qq=q(1)*(n-1:-1:1)/(n-1);
    else
        qq=(q(2):(-diff(q)/(n-2)):q(1)).*(n-1:-1:1)/(n-1);
    end

    [Wp,Wm]=MakeSMS(qq);


end

