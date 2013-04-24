function [ output ] = TestCounterEx( n,varnames )
%TESTCLIQUE Summary of this function goes here
%   Detailed explanation goes here


% [wp,wm,w]=CounterEx(n);
[wp,wm,w]=DiffJump(n);
[wp,wm ] = ModelOpt( wp,wm,1);
% q=ones(1,n-1);
% q(n/2+1)=1/n;
% q(n/2)=0;
% [wp,wm,w]=MakeSMS(q);
% [dwp,dwm]=MakePert(n,n/2,n);
% wp=wp+dwp;
% wm=wm+dwm;

% q=ones(1,n-1);
% q(2:(n/2-1))=1/2;
% [wp,wm,w]=MakeSMS(q);

% [wp,wm,w]=AllMtoP(n);

% if rcond(ones(n)-(wp+wm)/2) < 1e-15
%     crq=NaN;
%     crqm=NaN;
%     c2q=NaN;
%     Tmax=NaN;
%     return;
% end

[q,c]=SpectrumWpm(wp,wm,0.5,w);
% crqm=max(real(c.*sqrt(q)));
% crq=real(c'*sqrt(q));
% c2q=real((c.^2)'*q);
% init=q'*c;
% cq=max(real(q.*c));
% ev=q(2);
% 
% recqaqb=1./conj(q*ones(1,n)+ones(n,1)*q');
% recqaqb(1,:)=0;
% recqaqb(:,1)=0;
% L2=(c.*q)'*recqaqb*(c.*q);
% 
% T=FPT((wp+wm)/2);
% Tmax=max(max(T));
% 
% p=EqProb((wp+wm)/2);
% Eqbnd=p(n/2);
% 
% AI=real(sum(c)*c'*q);


output=zeros(length(varnames),1);
%
[tf,loc]=ismember('crqm',varnames);
if tf
    output(loc)=max(real(c.*sqrt(q)));
end
%
[tf,loc]=ismember('crq',varnames);
if tf
    output(loc)=real(c'*sqrt(q));
end
%
[tf,loc]=ismember('c2q',varnames);
if tf
    output(loc)=real((c.^2)'*q);
end
%
[tf,loc]=ismember('init',varnames);
if tf
    output(loc)=q'*c;
end
%
[tf,loc]=ismember('cqm',varnames);
if tf
    output(loc)=max(real(q.*c));
end
%
[tf,loc]=ismember('ev',varnames);
if tf
   output(loc)=q(2);
end
%
[tf,loc]=ismember('L2',varnames);
if tf
    recqaqb=1./conj(q*ones(1,n)+ones(n,1)*q');
    recqaqb(1,:)=0;
    recqaqb(:,1)=0;
    output(loc)=(c.*q)'*recqaqb*(c.*q);
end
%
[tf,loc]=ismember('A',varnames);
if tf
    output(loc)=real(sum(c));
end
%
[tf,loc]=ismember('AI',varnames);
if tf
    output(loc)=real(sum(c)*c'*q);
end
%
[tf,loc]=ismember('Eqbnd',varnames);
if tf
    p=EqProb((wp+wm)/2);
    output(loc)=p(n/2);
end
%
[tf,loc]=ismember('Tmax',varnames);
if tf
    T=FPT((wp+wm)/2);
    output(loc)=max(max(T));
end
%
[tf,loc]=ismember('qPi',varnames);
if tf
    output(loc)=-w'*(wp+wm)*w/(2*n);
end
%
[tf,loc]=ismember('etaPi',varnames);
if tf
    Zinv=ones(size(wp))-(wp+wm)/2;
    output(loc)=w'*(Zinv\w)/(2*n);
end
% allvars={'crqm','crq','c2q','init','cqm','ev','L2','AI','Eqbnd','Tmax'};
end

