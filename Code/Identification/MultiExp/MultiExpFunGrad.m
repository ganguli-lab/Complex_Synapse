function [ negloglike,grad ] = MultiExpFunGrad( params,data )
%[negloglike,grad]=MULTIEXPFUNGRAD(params,data)  likelihood of DATA under
%multi-exponential distribution with parameters PARAMS
%   params = column [log(q); c(1:end-1)]
%       c(end)=1-sum(c(1:end-1))
%   data = row
%   like = prod_data sum_{c,q} c*q*exp(-q*data)
%   negloglike = -log(like)/length(data)
%   grad(i) = d(negloglike)/d(params(i))

T=length(data);
n=ceil(length(params)/2);

%q_i = e^z_i
q=exp(params(1:n));

%c_n = 1 - sum c_i
cq=params(n+1:end);
cq=[cq;1-sum(cq)];
%cq=c.*q
cq=cq.*q;

qt=q*data;
expqt=exp(-qt);
onemat=ones(size(qt));

invlike = 1./(expqt'*cq);
negloglike=sum(log(invlike))/T;

gradc=diag(q)*expqt;
gradc=gradc-ones(n,1)*gradc(end,:);
gradc(end,:)=[];
gradz=diag(cq)*(onemat-qt).*expqt;
grad=-[gradz;gradc]*invlike/T;
% grad=sum(grad,2)/length(data);

end

