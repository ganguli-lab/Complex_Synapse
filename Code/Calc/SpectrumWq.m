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

Zinv=ones(n) - W;

[u,qa]=eig(-W);

% if rcond(u)<0.0001
%     disp('bad u');
% end
% if rcond(Zinv)<0.0001
%     disp('bad Z');
% end
ca = 2*fp*(1-fp) * (u\w) * sum((Zinv\q) * (Zinv\u), 1);


qa=diag(qa);
ca=diag(ca);

[qa,ix]=sort(qa);
ca=ca(ix);
end

