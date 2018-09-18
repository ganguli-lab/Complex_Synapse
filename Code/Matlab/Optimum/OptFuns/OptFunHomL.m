function [ f ] = OptFunHomL( x,s,fp,w,varargin )
%OPTFUNHOML function for fmincon on Laplace transform
%   x = parameters (off diagonal matrix elements)
%   f = function to be minimised (-SNR(t))
%   G = df/dx

[Wp,Wm,Q]=Params2MatsHom(x);

Wp=Wp + 1/(2*fp) * Q;
Wm=Wm + 1/(2*(1-fp)) * Q;



f=-SNRlaplace(s,Wp,Wm,fp,w);

% f=-real(crq(Wp,Wm,fp,w));
% f=-SNRarea( Wp, Wm, fp, w );
end

