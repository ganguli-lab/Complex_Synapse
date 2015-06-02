function [ A ] = OptFunChainL( qv,s,fp )
%OPTFUNCHAINL Laplace transform of SNR curve for asymmetric serial model for
%fmincon
%   Detailed explanation goes here

qp=qv(1:length(qv)/2);
qm=qv(length(qv)/2+1:end);
[Wp,Wm,w]=MakeMultistate(qp,qm);

p=fp/(1-fp)*qp./qm;
p=[1 cumprod(p)];
p=p/sum(p);

q=Wp-Wm;
A=2*fp*(1-fp)* p*q * ((s*eye(length(Wp))+ones(length(Wp))-fp*Wp-(1-fp)*Wm)\w);

A=-A;

end

