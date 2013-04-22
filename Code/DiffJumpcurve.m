function s=DiffJumpcurve(t,n)
error(CheckSize(t,@isrow))
error(CheckSize(n,@isscalar))
error(CheckValue(n,@(x)mod(x,2)==0,'even'))

[Wp,Wm,w]=DiffJump(n);
s=SNRcurve(t,Wp,Wm,0.5,w);















end