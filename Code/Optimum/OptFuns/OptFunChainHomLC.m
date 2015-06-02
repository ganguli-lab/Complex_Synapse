function [ A ] = OptFunChainHomLC( qv,s,fp )
%OPTFUNGRADCHAINHOMLC Laplace transform of SNR curve for asymmetric serial
%model with homeostatic plasticity for fmincon
%   Detailed explanation goes here

n=length(qv)/4;

q_pot_inc=qv(1:n);
q_pot_dec=qv(n+1:2*n);
q_dep_inc=qv(2*n+1:3*n);
q_dep_dec=qv(3*n+1:end);

Wp=StochastifyC(diag(q_pot_inc,1)+diag(q_pot_dec,-1));
Wm=StochastifyC(diag(q_dep_inc,1)+diag(q_dep_dec,-1));

w=BinaryWeights(n+1);

p=(fp*q_pot_inc+(1-fp)*q_dep_inc)./(fp*q_pot_dec+(1-fp)*q_dep_dec);
p=[1 cumprod(p)];
p=p/sum(p);

q=Wp-Wm;

Zinv=ones(length(Wp))-fp*Wp-(1-fp)*Wm;
Zinvs=s*eye(length(Wp))+Zinv;

a=q*(Zinvs\w);

A=2*fp*(1-fp)*p*a;

A=-A;

end

