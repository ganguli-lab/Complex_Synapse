function [ output ] = CounterExample( n )
%TESTCLIQUE Summary of this function goes here
%   Detailed explanation goes here

error(CheckSize(n,@isscalar));
error(CheckValue(n,@isint));
error(CheckValue(n,@(x) mod(x,2)==0,'even'));

fp=0.5;
qv=ones(1,n-1);
qv(n/2+1)=1/n;
% qv(n/2+1:end)=1/2;
[Wp,Wm,w]=MakeSMS(qv);

q=Wp-Wm;
W=Wm+fp*q;

ps=(w>0)';
ms=(w<0)';
bndp=find(ps,1,'first');
bndm=find(ms,1,'last');

p=EqProb(W);
T=FPT(W);

Tm=T(:,bndm);
Tp=T(:,bndp);

sig=p*q;

sigp2T=abs(sig(ps & sig>0)).*Tm(ps & sig>0)';
sigm2T=abs(sig(ms & sig<0)).*Tp(ms & sig<0)';

output=[sigm2T sigp2T];

% [qa,ca]=Spectrum(Wp,Wm,0.5,w);
% crq=max(real(ca.*sqrt(qa)));
% crq=ca'*sqrt(qa);
% init=qa'*ca;
% cq=max(real(qa.*ca));
% ev=qa(2);
% 
% recqaqb=1./conj(qa*ones(1,n)+ones(n,1)*qa');
% recqaqb(1,:)=0;
% recqaqb(:,1)=0;
% L2=(ca.*qa)'*recqaqb*(ca.*qa);
% 
% T=FPT((Wp+Wm)/2);
% Tmax=max(max(T));
% 
% p=EqProb((Wp+Wm)/2);
% Eqbnd=p(n/2);
% 
% AI=sum(ca)*ca'*qa;

end

