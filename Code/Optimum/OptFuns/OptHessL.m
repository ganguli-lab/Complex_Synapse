function [ hess ] = OptHessL( x,s,fp,w )
%[hess]=OPTHESSL(x,s,fp,w) Hessian of Laplace transform of SNR curve
%   x = parameters (off diagonal matrix elements)
%   f = function to be minimised (-SNR(t))
%   hess = d^2f/dx^2

[Wp,Wm] = Params2Mats(x);

[hesspp,hesspm,hessmp,hessmm] = SNRLaplaceHess(s,Wp,Wm,fp,w);

hess = - Mats2ParamsHess(hesspp,hesspm,hessmp,hessmm);

end

