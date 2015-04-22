function [ A ] = OptFunChainL( qv,s )
%OPTFUNCHAINL Laplace transform of SNR curve for asymmetric serial model for
%fmincon
%   Detailed explanation goes here

qp=qv(1:length(qv)/2);
qm=qv(length(qv)/2+1:end);
[Wp,Wm,w]=MakeMultistate(qp,qm);

p=qp./qm;
p=[1 cumprod(p)];
p=p/sum(p);

q=Wp-Wm;
A=0.5 * p*q * ((s*eye(length(Wp))+ones(length(Wp))-Wm-0.5*q)\w);

A=-A;

end

