function [ Wp,Wm,w ] = CounterEx( n,eps )
%[Wp,Wm,w]=COUNTEREX(n,eps) Counterexample generator
%   Detailed explanation goes here

error(CheckSize(n,@isscalar));
error(CheckValue(n,@(x)mod(x,2)==0,'even'));
existsAndDefault('eps',1/n);
error(CheckSize(eps,@isscalar));

q=ones(1,n-1);
q(n/2+1)=eps;
% q(1)=1/n;
[Wp,Wm,w]=MakeSMS(q);


end

