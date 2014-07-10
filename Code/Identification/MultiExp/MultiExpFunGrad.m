function [ negloglike,grad ] = MultiExpFunGrad( params,data )
%[negloglike,grad]=MULTIEXPFUNGRAD(params,data)  likelihood of DATA under
%multi-exponential distribution with parameters PARAMS
%   params = column [q; c(1:end-1)]
%       c(end)=1-sum(c(1:end-1))
%   data = row
%   like = prod_data sum_{c,q} c*q*exp(-q*data)
%   negloglike = -log(like)
%   grad(i) = d(negloglike)/d(params(i))

n=ceil(length(params)/2);
q=params(1:n);
c=params(n+1:end);
c=[c;1-sum(c)];

qt=q*data;
expqt=exp(-qt);
onemat=ones(size(qt));

like = (c.*q)'*expqt;
negloglike=-sum(log(like));

gradc=diag(q)*expqt;
gradc=gradc-ones(n,1)*gradc(end,:);
gradc(end,:)=[];
gradq=diag(c)*(onemat-qt).*expqt;
grad=-[gradq;gradc]*diag(1./like);
grad=sum(grad,2);

end

