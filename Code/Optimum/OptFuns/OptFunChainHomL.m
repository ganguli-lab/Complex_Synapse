function [ A ] = OptFunChainHomL( qv,s )
%OPTFUNGRADCHAINHOML Laplace transform of SNR curve for asymmetric serial
%model with homeostatic plasticity for fmincon
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

p=(qp+qhp)./(qm+qhm);
p=[1 cumprod(p)];
p=p/sum(p);

q=Wp-Wm;

Zinv=ones(length(Wp))-Wm-0.5*q-H;
Zinvs=s*eye(length(Wp))+Zinv;

a=q*(Zinvs\w);

A=0.5*p*a;

A=-A;

end

