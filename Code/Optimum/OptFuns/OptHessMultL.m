 function [ hessV ] = OptHessMultL( x0,s,fp,w,xp )
%[hess]=OPTHESSL(x0,s,fp,w,xp) multiplies diff of transtion matrices by
%hessian of Laplace transform of SNR curve
%   x0 = parameters (off diagonal matrix elements) where hessian evaluated
%   xp = parameters (off diagonal matrix elements) to multiply by hessian
%   f = function to be minimised (-SNR(t))
%   hess = d^2f/dx^2

[Wp,Wm] = Params2Mats(x0);

[Vp,Vm] = Params2Mats(xp);

fm=1-fp;

ev=ones(1,length(Wp));
q=Wp-Wm;
W=Wm+fp*q;

Zinv=ones(length(Wp))-W;

p=ev/Zinv;

Zinvs=s*eye(length(Wp))+Zinv;

Zw=Zinvs\w;
a=q*Zw;
Za=Zinv\a;
c=(p*q)/Zinvs;

Zap=Za*p;
Zwp=Zw*p;
Zwc=Zw*c;

ZqZ = (Zinv \ q) / Zinvs;

hess2p  = Zinv  \ Vp * Zap + Zap * Vp / Zinv;
hess4p  = ZqZ   * Vp * Zwp + Zwp * Vp * ZqZ;
hess10p = Zinvs \ Vp * Zwc + Zwc * Vp / Zinvs;
hess3p  = Zinv  \ Vp * Zwp + Zwp * Vp / Zinvs;
hess6p  = Zinvs \ Vp * Zwp + Zwp * Vp / Zinv;

hess2m  = Zinv  \ Vm * Zap + Zap * Vm / Zinv;
hess4m  = ZqZ   * Vm * Zwp + Zwp * Vm * ZqZ;
hess10m = Zinvs \ Vm * Zwc + Zwc * Vm / Zinvs;
hess3m  = Zinv  \ Vm * Zwp + Zwp * Vm / Zinvs;
hess6m  = Zinvs \ Vm * Zwp + Zwp * Vm / Zinv;

hessVpp = fp^2  * (hess2p + hess4p + hess10p) + fp * (hess3p +      hess6p);
hessVmp = fp*fm * (hess2p + hess4p + hess10p) + fm *  hess3p - fp * hess6p;
hessVpm = fp*fm * (hess2m + hess4m + hess10m) - fp *  hess3m + fm * hess6m;
hessVmm = fm^2  * (hess2m + hess4m + hess10m) - fp * (hess3m +      hess6m);

hessVp = (hessVpp + hessVpm)';
hessVm = (hessVmp + hessVmm)';

hessVp = hessVp - diag(hessVp) * ev;
hessVm = hessVm - diag(hessVm) * ev;

hessV = - Mats2Params(hessVp,hessVm);

end

