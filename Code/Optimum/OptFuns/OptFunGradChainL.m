function [ A,gr ] = OptFunGradChainL( qv,s )
%OPTFUNGRADCHAINL Laplace transform of SNR curve for symmetric serial model
%and gradient for fmincon
%   Detailed explanation goes here

[Wp,Wm,w]=MakeSMS(qv);

p=qv./wrev(qv);
p=[1 cumprod(p)];
p=p/sum(p);

q=Wp-Wm;

Zinv=ones(length(Wp))-Wm-0.5*q;
Zinvs=s*eye(length(Wp))+Zinv;

a=0.5*q*(Zinvs\w);
c=(p*q)/Zinvs;

A=0.5*p*a;

% %dA(s)/dq_(i,i+1)
% dAdq = p(1:end-1).*(diff(Zinvs\w))';
% %dA(s)/dW^F_(i,i+1)
% dAdW = p(1:end-1) .* diff(Zinv\a)' + c(1:end-1).*diff(Zinvs\w)';
% %dA/dWp_(i,i+1)+dA/dWm_(n-i+1,n-i)
% gr = 2*dAdq+dAdW;

%dA/dWp_(i,i+1)+dA/dWm_(n-i+1,n-i)
gr= p(1:end-1).*diff(Zinv\a)' + (p(1:end-1)+0.5*c(1:end-1)).*diff(Zinvs\w)';

A=-A;
gr=-gr;

end

