function [ crq,init,cq,ev,Tmax,Eqbnd,L2,AI ] = TestClique( n )
%TESTCLIQUE Summary of this function goes here
%   Detailed explanation goes here

assert(isscalar(n));
assert(mod(n,6)==0);

% [wp,wm,w]=CliqueChain(n);

q=ones(1,n-1);
q(n/2+1)=1/n;
[wp,wm,w]=MakeSMS(q);

% [wp,wm,w]=AllMtoP(n);


[q,c]=Spectrum(wp,wm,0.5,w);
% crq=max(real(c.*sqrt(q)));
crq=c'*sqrt(q);
init=q'*c;
cq=max(real(q.*c));
ev=q(2);

recqaqb=1./conj(q*ones(1,n)+ones(n,1)*q');
recqaqb(1,:)=0;
recqaqb(:,1)=0;
L2=(c.*q)'*recqaqb*(c.*q);

T=FPT((wp+wm)/2);
Tmax=max(max(T));

p=EqProb((wp+wm)/2);
Eqbnd=p(n/2);

AI=sum(c)*c'*q;

end

