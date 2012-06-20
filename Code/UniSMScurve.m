function s=UniSMScurve(t,n,q,eps)
existsAndDefault('eps',q);
error(CheckSize(t,@isrow))
error(CheckSize(n,@isscalar))
error(CheckSize(q,@isscalar))
error(CheckSize(eps,@isscalar))
error(CheckValue(n,@(x)mod(x,2)==0,'even'))
error(CheckValue(q,@(x)inrange(x,0,1),'inrange(0,1)'))
error(CheckValue(eps,@(x)inrange(x,0,1),'inrange(0,1)'))

qq=q*ones(1,n-1);
qq(1)=eps;
[Wp,Wm,w]=MakeSMS(qq);
s=SNRcurve(t,Wp,Wm,0.5,w);















end