function [ negloglike ] = MultiExpFun( params,data )
%negloglike=MULTIEXPFUN(PARAMS,DATA) likelihood of DATA under
%multi-exponential distribution with parameters PARAMS
%   params = column [log(q); c(1:end-1)]
%       c(end)=1-sum(c(1:end-1))
%   data = row
%   like = prod_data sum_{c,q} c*q*exp(-q*data)
%   negloglike = -log(like)/length(data)

n=ceil(length(params)/2);
q=exp(params(1:n));
c=params(n+1:end);
c=[c;1-sum(c)];

like = (c.*q)'*exp(-q*data);
negloglike=-sum(log(like))/length(data);


end

