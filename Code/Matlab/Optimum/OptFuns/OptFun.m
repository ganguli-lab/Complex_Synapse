function [ f ] = OptFun( x,t,fp,w,varargin )
%OPTFUN function and gradient for fmincon
%   x = parameters (off diagonal matrix elements)
%   f = function to be minimised (-SNR(t))
%   G = df/dx

[Wp,Wm]=Params2Mats(x);

q=Wp-Wm;
W=fp*q + Wm;
p=EqProb(W,varargin{:});

f=-p*q*expm(W*t)*w;

% f=-real(crq(Wp,Wm,fp,w));
% f=-SNRarea( Wp, Wm, fp, w );
end

