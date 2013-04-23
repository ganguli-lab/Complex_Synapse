function [ qa,ca ] = SpectrumWq( W,q,fp,w,varargin )
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
% assert(mod(n,2)==0)

DegThresh=1e-3;
RCondThresh=1e-5;
varargin=assignApplicable(varargin);

% w=ones(n,1);
% w(1:(n/2))=-1;


[u,qa]=eig(-W);
qa=diag(qa);
[qa,ix]=sort(qa);
u=u(:,ix);

%this method doesn't work when eigenvalues are nearly degenerate
if min(diff(qa)) > DegThresh
    [v,qb]=eig(-W');
    qb=diag(qb);
    [~,ix]=sort(qb);
    v=v(:,ix);
    v=diag(1./diag(v.'*u)) * v.';

    Z=u * diag(1./[1;qa(2:end)]) * v;
    p=v(1,:);
    p=p/sum(p);
    ca = 2*fp*(1-fp) * (v*w) .* (p*q*Z*u)';
    return;
end


%this method doesn't work when eigenvectors are nearly parallel
% if rcond(u) > RCondThresh
    Zinv=ones(n) - W;
    %this method doesn't work when eigenvectors are nearly parallel or W
    %non-ergodic
    if rcond(Zinv) > RCondThresh
        ca = 2*fp*(1-fp) * (u\w) * sum((Zinv\q) * (Zinv\u), 1);
        ca=diag(ca);
    else
        Z=u * diag(1./[1;qa(2:end)]) / u;
        p=[1 zeros(1,length(qa)-1)] / u;
        p=p/sum(p);
        ca = 2*fp*(1-fp) * (u\w) .* (p*q*Z*u)';
    end
% else
%     qa=NaN(n,1);
%     ca=qa;
% end

end

