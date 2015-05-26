function [ A,gr ] = OptFunGradChainHomL( qv,s )
%OPTFUNGRADCHAINHOML Laplace transform of SNR curve and gradient for
%asymmetric serial model with homeostatic plasticity for fmincon
%   Detailed explanation goes here

n=length(qv)/4;

%activity dependent plasticity
qp=qv(1:n);
qm=qv(n+1:2*n);
[Wp,Wm,w]=MakeMultistate(qp,qm);
%acticity independent plasticity
qhp=qv(2*n+1:3*n);
qhm=qv(3*n+1:end);
[Hp,Hm]=MakeMultistate(qhp,qhm);
H=Hp+Hm;

p=(0.5*qp+qhp)./(0.5*qm+qhm);
p=[1 cumprod(p)];
p=p/sum(p);

q=Wp-Wm;

Zinv=ones(length(Wp))-Wm-0.5*q-H;
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

%dA/dWp_(i,i+1)+dA/dWm_(i,i+1)
grhp= 0.5*p(1:end-1).*Za + 0.5*c(1:end-1).*Zw;
%dA/dWp_(i+1,i)+dA/dWm_(i+1,i)
grhm= -0.5*p(2:end).*Za -0.5*c(2:end).*Zw;
%dA/dWp_(i,i+1)
grp= grhp + p(1:end-1).*Zw;
%dA/dWm_(i+1,i)
grm= grhm + p(2:end).*Zw;

A=-A;
gr=-0.5*[grp grm grhp grhm];

end

