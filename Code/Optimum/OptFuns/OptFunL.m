function [ f ] = OptFunL( x,s,fp,w,varargin )
%OPTFUNL function and gradient for fmincon on Laplace transform
%   x = parameters (off diagonal matrix elements)
%   f = function to be minimised (-SNR(t))
%   G = df/dx

[Wp,Wm]=Params2Mats(x);


f=-SNRlaplace(s,Wp,Wm,fp,w);

% f=-real(crq(Wp,Wm,fp,w));
% f=-SNRarea( Wp, Wm, fp, w );
end

