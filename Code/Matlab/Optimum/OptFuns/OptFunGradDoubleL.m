function [ f,gr ] = OptFunGradDoubleL( x,s,fp,w,varargin )
%OPTFUNGRADDOUBLEL function and gradient for fmincon on Laplace transform
%   x = parameters (off diagonal matrix elements)
%   f = function to be minimised (-SNR(t))
%   G = df/dx

[Wp,Wm]=Params2Mats(x);

[ A,dAdWp,dAdWm ] = DoubleLaplace( s,Wp,Wm,fp,w );

f=-A;
gr = -Mats2Params(dAdWp,dAdWm);


end

