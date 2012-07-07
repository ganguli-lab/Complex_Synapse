function [ qa,ca ] = SpectrumWq( W,q,fp, w )
%[QA,CA]=SPECTRUM(W,q,w) Weights and eigenvalues of SNR curve for complex synapse (cts time)
%   W = forgetting transition rates: f^+W^+ + f^-W^-
%   q = encoding transition rates: W^+ - W^-
%   w = Weights of states (+/-1)
%   QA = eigenvalues (decay rate)
%   CA = weights (contribution to area)

error(CheckSize(W,@ismat));%matrix
error(CheckSize(W,@issquare));%square
error(CheckSize(q,@(x)samesize(W,x),'samesize(W)'));%also square matrix of same size
error(CheckSize(w,@iscol));

n=size(W,1);
assert(mod(n,2)==0)

% w=ones(n,1);
% w(1:(n/2))=-1;


[u,qa]=eig(-W);
qa=diag(qa);
[qa,ix]=sort(qa);
u=u(:,ix);

[v,qb]=eig(-W');
qb=diag(qb);
[~,ix]=sort(qb);
v=conj(v(:,ix));
v=diag(1./diag(v'*u)) * v';

% Zinv=ones(n) - W;
% ca = 2*fp*(1-fp) * (u\w) * sum((Zinv\q) * (Zinv\u), 1);
% ca=diag(ca);
% ca=ca(ix);

% if rcond(u) < 1e-7
%     qa=NaN(n,1);
%     ca=qa;
%     return;
% end

% Z=u * diag(1./[1;qa(2:end)]) / u;
% p=[1 zeros(1,length(qa)-1)] / u;
% p=p/sum(p);
% ca = 2*fp*(1-fp) * (u\w) .* (p*q*Z*u)';

Z=u * diag(1./[1;qa(2:end)]) * v;
p=v(1,:);
p=p/sum(p);
ca = 2*fp*(1-fp) * (v*w) .* (p*q*Z*u)';

end

