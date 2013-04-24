function [ f ] = OptFun( x,t,fp,w )
%OPTFUN function and gradient for fmincon
%   x = parameters (off diagonal matrix elements)
%   f = function to be minimised (-SNR(t))
%   G = df/dx

[Wp,Wm]=Params2Mats(x);

% q=Wp-Wm;
% W=fp*q + Wm;
% p=EqProb(W);
% 
% f=-p*q*expm(W*t)*w;

f=-real(crq(Wp,Wm,fp,w));

end

