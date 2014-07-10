function [ negloglike,grad,hess ] = MultiExpFunGradHess( params,data )
%[negloglike,grad,hess]=MULTIEXPFUNGRADHESS(params,data)  likelihood of
%DATA under multi-exponential distribution with parameters PARAMS
%   params = column [q; c(1:end-1)]
%       c(end)=1-sum(c(1:end-1))
%   data = row
%   like = prod_data sum_{c,q} c*q*exp(-q*data)
%   negloglike = -log(like)
%   grad(i) = d(negloglike)/d(params(i))
%   hess(i,j) = d^2(negloglike)/d(params(i))d(params(j))

n=ceil(length(params)/2);
q=params(1:n);
c=params(n+1:end);
c=[c;1-sum(c)];
T=length(data);

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

hesscc=mmx('mult', reshape(gradc,n-1,1,T), reshape(gradc,1,n-1,T) );

hesscq = mmx('mult', reshape(gradc,n-1,1,T), reshape(gradq,1,n,T) );
%needs fixing

hessqq = mmx('mult', reshape(gradq,n,1,T), reshape(gradq,1,n,T) );
%needs fixing

hess = cat(1, cat(2, hesscc, hesscq ), cat(2, permute(hesscq,[2 1 3]), hessqq ) );
end

