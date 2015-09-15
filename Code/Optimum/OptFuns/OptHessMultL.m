function [ hessV ] = OptHessMultL( x0,s,fp,w,xp )
%[hess]=OPTHESSL(x0,s,fp,w,xp) multiplies diff of transtion matrices by hessian
%of Laplace transform of SNR curve
%   x0 = parameters (off diagonal matrix elements)
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

hess2p  = Zinv  \ Vp * Zap;
hess4p  = ZqZ   * Vp * Zwp;
hess10p = Zinvs \ Vp * Zwc;
hess3p  = Zinv  \ Vp * Zwp;
hess6p  = Zinvs \ Vp * Zwp;

hess2m  = Zinv  \ Vm * Zap;
hess4m  = ZqZ   * Vm * Zwp;
hess10m = Zinvs \ Vm * Zwc;
hess3m  = Zinv  \ Vm * Zwp;
hess6m  = Zinvs \ Vm * Zwp;

hessVpp = fp^2  * (hess2p + hess4p * hess10p) + fp * (hess3p + hess6p);
hessVmp = fp*fm * (hess2p + hess4p * hess10p) + fm * hess3p - fp * hess6p;
hessVpm = fp*fm * (hess2m + hess4m * hess10m) - fp * hess3m + fm * hess6m;
hessVmm = fm^2  * (hess2m + hess4m * hess10m) - fp * (hess3m + hess6m);

hessVp = hessVpp + hessVpm;
hessVm = hessVmp + hessVmm;

hessV = - Mats2Params(hessVp,hessVm);

end

