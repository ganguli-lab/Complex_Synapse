function [ A,gr ] = OptFunGradChainHomL( qv,s,fp )
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

p=(fp*qp+qhp)./((1-fp)*qm+qhm);
p=[1 cumprod(p)];
p=p/sum(p);

Wp=Wp+0.5/fp*H;
Wm=Wm+0.5/(1-fp)*H;

q=Wp-Wm;

Zinv=ones(length(Wp))-fp*Wp-(1-fp)*Wm;
Zinvs=s*eye(length(Wp))+Zinv;

a=q*(Zinvs\w);
c=(p*q)/Zinvs;

Za=diff(Zinv\a)';
Zw=diff(Zinvs\w)';

A=2*fp*(1-fp)*p*a;

% %dA(s)/dq_(i,i+1)
% dAdq = p(1:end-1).*(diff(Zinvs\w))';
% %dA(s)/dW^F_(i,i+1)
% dAdW = p(1:end-1) .* diff(Zinv\a)' + c(1:end-1).*diff(Zinvs\w)';
% %dA/dWp_(i,i+1)+dA/dWm_(n-i+1,n-i)
% gr = 2*dAdq+dAdW;

%dA/dWp_(i,i+1)*dA/dWm_(i,i+1)
grfp= p(1:end-1).*Za + c(1:end-1).*Zw;
%dA/dWp_(i+1,i)*dA/dWm_(i+1,i)
grfm= -p(2:end).*Za -c(2:end).*Zw;
%fm*dA/dWp_(i,i+1)-fp*dA/dWm_(i,i+1)
grqp=  p(1:end-1).*Zw;
%fm*dA/dWp_(i+1,i)-fp*dA/dWm_(i+1,i)
grqm= p(2:end).*Zw;

grp=fp*grfp+grqp;
grm=(1-fp)*grfm-grqm;
grhp=grp+0.5*(1/fp-1/(1-fp))*grqp;
grhm=grm+0.5*(1/fp-1/(1-fp))*grqm;

A=-A;
gr=-2*fp*(1-fp)*[grp grm grhp grhm];

end

