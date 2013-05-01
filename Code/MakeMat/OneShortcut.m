function [ Wp,Wm,w ] = OneShortcut( n,i,j,qs,qd )
%[Wp,Wm,w]=ONESHORTCUT(n,i,j,qs,qd) Synapse with one shortcut for
%potentiation, flip symmetric for depression
%   n  = # states
%   i  = starting state for shortcut
%   j  = final state for shortcut
%   qs = transition prob of shortcut (default: 1)
%   qd = transition prob of direct link out of state i (default: 1-qs)

error(CheckSize(n,@isscalar));
error(CheckSize(i,@isscalar));
error(CheckSize(j,@isscalar));
error(CheckValue(n,@isint));
error(CheckValue(i,@isint));
error(CheckValue(j,@isint));
error(CheckValue(i,@(x)inrange(x,1,n),'inrange(1,n)'));
error(CheckValue(j,@(x)inrange(x,1,n),'inrange(1,n)'));

existsAndDefault('qs',1);
error(CheckSize(qs,@isscalar));
error(CheckValue(qs,@(x)inrange(x,0,1),'inrange(0,1)'));
existsAndDefault('qd',1-qs);
error(CheckSize(qd,@isscalar));
error(CheckValue(qd,@(x)inrange(x,0,1-qs),'inrange(0,1-qs)'));

q=ones(1,n-1);
q(i)=qd;

[Wp,Wm,w]=MakeSMS(q);

Wp(i,j)=qs;
Wm(n-i+1,n-j+1)=qs;
Wp=StochastifyC(Wp);
Wm=StochastifyC(Wm);

end

