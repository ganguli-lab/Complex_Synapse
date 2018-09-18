function [ genfun ] = MultiExpGenFunTrue( s,params )
%genfun=MULTIEXPGENFUN(s,data) estimate of moment generating function
%   genfun = <exp(-st)> = sum c*q/(s+q)
%   s      = column of s values
%   params = column [log(q); c(1:end-1)]
%       c(end)=1-sum(c(1:end-1))

n=ceil(length(params)/2);
q=exp(params(1:n));
c=params(n+1:end);
c=[c;1-sum(c)];

den=1./(bsxfun(@plus,s,q'));
genfun=den*(c.*q);


end

