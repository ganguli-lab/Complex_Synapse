function [ A ] = SNRrunAve( t,Wp,Wm,fp,w )
%A=SNRLAPLACE(s,Wp,Wm,fp,w) running maverage of SNR curve for complex
%synapse (cts time), i.e. (Laplace Transform @ s=1/t)/t.
%   A(s) = int exp(-s*t)*SNR(t) dt
%   t  = time = 1/parameter of Laplace transform
%   WP = potentiation transition rates
%   WM = depression transition rates
%   FP = Fraction of potentiation transitions
%   w  = Weights of states (+/-1)
error(CheckSize(Wp,@ismat));%matrix
error(CheckSize(Wp,@issquare));%square
error(CheckSize(Wm,@samesize,'samesize(Wp)',Wp));%also square matrix of same size
error(CheckSize(fp,@isscalar));
error(CheckValue(fp,@inrange,'inrange(0,1)',0,1));%fp in [0,1]
error(CheckSize(w,@iscol));
error(CheckValue(w,@(x) all(x.^2==1),'all w = +/-1'));
error(CheckSize(w,@samelength,'samelength(Wp)',Wp));%same size

A = SNRlaplace(1./t,Wp,Wm,fp,w)./t;


end

