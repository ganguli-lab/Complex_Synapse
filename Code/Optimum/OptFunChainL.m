function [ A ] = OptFunChainL( qv,s )
%OPTFUNCHAINL Laplace transform of SNR curve for symmetric serial model for
%fmincon
%   Detailed explanation goes here

[Wp,Wm,w]=MakeSMS(qv);

p=qv./wrev(qv);
p=[1 cumprod(p)];
p=p/sum(p);

q=Wp-Wm;
A=0.5 * p*q * ((s*eye(length(Wp))+ones(length(Wp))-Wm-0.5*q)\w);

A=-A;

end

