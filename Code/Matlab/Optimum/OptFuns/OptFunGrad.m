function [ f,G ] = OptFunGrad( x,t,fp,w,varargin )
%[f,G]=OPTFUNGRAD(x) function and gradient for fmincon
%   x = parameters (off diagonal matrix elements)
%   f = function to be minimised (-SNR(t))
%   G = df/dx

[Wp,Wm]=Params2Mats(x);

[f,dSdWp,dSdWm] = GradSNR(t,Wp,Wm,fp,w,varargin{:});
f=-f;
G=Mats2Params(dSdWp,dSdWm);
G=-G;


end

