function [Wp,Wm,w]= CascadeOriginal( xp,xm,n,lambda )
%[Wp,Wm,w]=CASCADEORIGINAL(xp,xm,n) Makes Cascade model, as originally
%designed
%   XP,XM  = ratio of nearby transition rates
%   n  = number of states (even)
%   WP = potentiation transition rates
%   WM = depression transition rates
%   w  = Weights of states (+/-1)

error(CheckSize(xp,@isscalar));
error(CheckValue(xp,@(x) inrange(x,0,0.5)));
error(CheckSize(xm,@isscalar));
error(CheckValue(xm,@(x) inrange(x,0,0.5)));
error(CheckSize(n,@isscalar));
error(CheckValue(n,@isint));
existsAndDefault('lambda',1);

qp=(xp.^abs((1:n-1)-n/2))/(1-xp);
qm=(xm.^abs((1:n-1)-n/2))/(1-xm);
[Wp,~,w]=CascadeMSinterp(qp,qp,0.5,lambda);
[~,Wm]=CascadeMSinterp(qm,qm,0.5,lambda);

end

