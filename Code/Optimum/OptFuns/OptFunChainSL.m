function [ A ] = OptFunChainSL( qv,s,fp )
%OPTFUNCHAINL Laplace transform of SNR curve for symmetric serial model for
%fmincon
%   Detailed explanation goes here

[Wp,Wm,w]=MakeSMS(qv);

p=fp*qv./wrev((1-fp)*qv);
p=[1 cumprod(p)];
p=p/sum(p);

q=Wp-Wm;
A=2*fp*(1-fp)* p*q * ((s*eye(length(Wp))+ones(length(Wp))--fp*Wp-(1-fp)*Wm)\w);

A=-A;

end

