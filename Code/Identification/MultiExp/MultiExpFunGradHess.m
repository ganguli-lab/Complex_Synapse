function [ negloglike,grad,hess ] = MultiExpFunGradHess( params,data )
%[negloglike,grad,hess]=MULTIEXPFUNGRADHESS(params,data)  likelihood of
%DATA under multi-exponential distribution with parameters PARAMS
%   params = column [q; c(1:end-1)]
%       c(end)=1-sum(c(1:end-1))
%   data = row
%   like = prod_data sum_{c,q} c*q*exp(-q*data)
%   negloglike = -log(like)
%   grad(i) = d(negloglike)/d(params(i))/length(data)
%   hess(i,j) = d^2(negloglike)/d(params(i))d(params(j))

T=length(data);
n=ceil(length(params)/2);

%q_i = e^z_i
q=exp(params(1:n));

%c_n = 1 - sum c_i
cq=params(n+1:end);
cq=[cq;1-sum(cq)];
%cq=c.*q
cq=cq.*q;

% %c_n = 1 - sum c_i
% c=params(n+1:end);
% c=[c;1-sum(c)];
% %cq=c.*q
% cq=c.*q;

qt=q*data;
expqt=exp(-qt);
onemat=ones(size(qt));
onev=ones(n,1);

invlike = 1./(expqt'*cq);
negloglike=sum(log(invlike))/T;

gradc=diag(q)*expqt;
gradc=gradc-ones(n,1)*gradc(end,:);
% gradc=gradc-onev*gradc(end,:);
gradc(end,:)=[];
gradz=diag(cq)*(onemat-qt).*expqt;
grad=-[gradz;gradc]*invlike/T;
% grad=sum(grad,2)/T;

% hesscc=mmx('mult', reshape(gradc,n-1,1,T), reshape(gradc,1,n-1,T) );
hesscc=gradc*diag(invlike.^2)*gradc';

% hesscq = mmx('mult', reshape(gradc,n-1,1,T), reshape(gradq,1,n,T) );
hesscz=gradc*diag(invlike.^2)*gradz';
%needs fixing
% hessczd = diag(1./c)*gradz*invlike;
% hessczd = hessczd - onev*hessczd(end);
% hessczd(end)=[];
% hesscz=hesscz-diag(hessczd);

% hessqq = mmx('mult', reshape(gradz,n,1,T), reshape(gradz,1,n,T) );
hesszz=gradz*diag(invlike.^2)*gradz';
%needs fixing
% hesszz = hesszz  - diag(((diag(c)*qt).*(2*onemat-qt).*expqt)*invlike);

% hess = sum(cat(1, cat(2, hesscc, hesscz ), cat(2, permute(hesscz,[2 1 3]), hesszz ) ),3)/T;
hess = [hesszz hesscz'; hesscz hesscc]/T;

end

