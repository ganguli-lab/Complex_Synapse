function [ A,gr ] = OptFunGradChainL( qv,s )
%OPTFUNGRADCHAINL Laplace transform of SNR curve for asymmetric serial model
%and gradient for fmincon
%   Detailed explanation goes here

qp=qv(1:length(qv)/2);
qm=qv(length(qv)/2+1:end);
[Wp,Wm,w]=MakeMultistate(qp,qm);

p=qp./qm;
p=[1 cumprod(p)];
p=p/sum(p);

q=Wp-Wm;

Zinv=ones(length(Wp))-Wm-0.5*q;
Zinvs=s*eye(length(Wp))+Zinv;

a=q*(Zinvs\w);
c=(p*q)/Zinvs;

Za=diff(Zinv\a)';
Zw=diff(Zinvs\w)';

A=0.5*p*a;

% %dA(s)/dq_(i,i+1)
% dAdq = p(1:end-1).*(diff(Zinvs\w))';
% %dA(s)/dW^F_(i,i+1)
% dAdW = p(1:end-1) .* diff(Zinv\a)' + c(1:end-1).*diff(Zinvs\w)';
% %dA/dWp_(i,i+1)+dA/dWm_(n-i+1,n-i)
% gr = 2*dAdq+dAdW;

%dA/dWp_(i,i+1)+dA/dWm_(n-i+1,n-i)
grp= 0.5*p(1:end-1).*Za + (p(1:end-1)+0.5*c(1:end-1)).*Zw;
grm= -0.5*p(2:end).*Za + (p(2:end)-0.5*c(2:end)).*Zw;

A=-A;
gr=-0.5*[grp grm];

end

