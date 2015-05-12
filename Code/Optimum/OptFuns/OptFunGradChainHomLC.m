function [ A,gr ] = OptFunGradChainHomLC( qv,s )
%OPTFUNGRADCHAINHOML Laplace transform of SNR curve and gradient for
%asymmetric serial model with homeostatic plasticity for fmincon
%   Detailed explanation goes here

n=length(qv)/4;

q_pot_inc=qv(1:n);
q_pot_dec=qv(n+1:2*n);
q_dep_inc=qv(2*n+1:3*n);
q_dep_dec=qv(3*n+1:end);

Wp=StochastifyC(diag(q_pot_inc,1)+diag(q_pot_dec,-1));
Wm=StochastifyC(diag(q_dep_inc,1)+diag(q_dep_dec,-1));

w=BinaryWeights(n+1);

p=(q_pot_inc+q_dep_inc)./(q_pot_dec+q_dep_dec);
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

%dA/dW(i,i+1)
dAdW_inc= p(1:end-1).*Za + c(1:end-1).*Zw;
%dA/dW(i+1,i)
dAdW_dec= -p(2:end).*Za - c(2:end).*Zw;
%dA/dq_(i,i+1)
dAdq_inc= p(1:end-1).*Zw;
%dA/dq_(i+1,i)
dAdq_dec= - p(2:end).*Zw;

gr_pot_inc=0.5*dAdW_inc+dAdq_inc;
gr_pot_dec=0.5*dAdW_dec+dAdq_dec;
gr_dep_inc=0.5*dAdW_inc-dAdq_inc;
gr_dep_dec=0.5*dAdW_dec-dAdq_dec;

A=-A;
gr=-0.5*[gr_pot_inc gr_pot_dec gr_dep_inc gr_dep_dec];

end

