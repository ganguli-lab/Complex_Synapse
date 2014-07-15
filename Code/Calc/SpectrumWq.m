function [ qa,ca ] = SpectrumWq( W,q,fp,w,varargin )
%[QA,CA]=SPECTRUM(W,q,w) Weights and eigenvalues of SNR curve for complex synapse (cts time)
%   W = forgetting transition rates: f^+W^+ + f^-W^-
%   q = encoding transition rates: W^+ - W^-
%   w = Weights of states (+/-1)
%   QA = eigenvalues (decay rate)
%   CA = weights (contribution to area)

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='SpectrumWq';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addRequired('W',@(x)validateattributes(x,{'numeric'},{'2d','square'},'SpectrumWq','W',1));
    p.addRequired('q',@(x)validateattributes(x,{'numeric'},{'2d','square'},'SpectrumWq','q',2));
    p.addRequired('fp',@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1},'SpectrumWq','fp',3));
    p.addRequired('w',@(x)validateattributes(x,{'numeric'},{'column'},'SpectrumWq','w',4));
    p.addParameter('DegThresh',1e-3,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'},'SpectrumWq','DegThresh'));
    p.addParameter('RCondThresh',1e-5,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'},'SpectrumWq','RCondThresh'));
end
p.parse(W,q,fp,w,varargin{:});
r=p.Results;

error(CheckSize(r.q,@(x)samesize(r.W,x),'samesize(W)'));%also square matrix of same size
error(CheckSize(r.w,@(x)length(x)==length(r.W),'samesize(W)'));

n=size(r.W,1);
% assert(mod(n,2)==0)


% w=ones(n,1);
% w(1:(n/2))=-1;


[u,qa]=eig(-r.W);
qa=diag(qa);
[qa,ix]=sort(qa);
u=u(:,ix);

%this method doesn't work when eigenvalues are nearly degenerate
if min(diff(qa)) > r.DegThresh
    [v,qb]=eig(-r.W');
    qb=diag(qb);
    [~,ix]=sort(qb);
    v=v(:,ix);
    v=diag(1./diag(v.'*u)) * v.';

    Z=u * diag(1./[1;qa(2:end)]) * v;
    p=v(1,:);
    p=p/sum(p);
    ca = 2*r.fp*(1-r.fp) * (v*r.w) .* (p*r.q*Z*u)';
    return;
end


%this method doesn't work when eigenvectors are nearly parallel
% if rcond(u) > RCondThresh
    Zinv=ones(n) - r.W;
    %this method doesn't work when eigenvectors are nearly parallel or W
    %non-ergodic
    if rcond(Zinv) > r.RCondThresh
        ca = 2*fp*(1-fp) * (u\r.w) * sum((Zinv\r.q) * (Zinv\u), 1);
        ca=diag(ca);
    else
        Z=u * diag(1./[1;qa(2:end)]) / u;
        p=[1 zeros(1,length(qa)-1)] / u;
        p=p/sum(p);
        ca = 2*r.fp*(1-r.fp) * (u\r.w) .* (p*r.q*Z*u)';
    end
% else
%     qa=NaN(n,1);
%     ca=qa;
% end

end

