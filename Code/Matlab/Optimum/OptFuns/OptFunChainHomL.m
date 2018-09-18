function [ A ] = OptFunChainHomL( qv,s,fp )
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

p=(fp*qp+qhp)./((1-fp)*qm+qhm);
p=[1 cumprod(p)];
p=p/sum(p);

Wp=Wp+0.5/fp*H;
Wm=Wm+0.5/(1-fp)*H;

q=Wp-Wm;

Zinv=ones(length(Wp))-fp*Wp-(1-fp)*Wm;
Zinvs=s*eye(length(Wp))+Zinv;

a=q*(Zinvs\w);

A=2*fp*(1-fp)*p*a;

A=-A;

end

