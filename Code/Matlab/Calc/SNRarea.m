function [ area ] = SNRarea( Wp, Wm, fp, w )
%AREA=SNRAREA(WP,WM,FP,w) Area under SNR curve for complex synapse (cts time)
%   WP = potentiation transition rates
%   WM = depression transition rates
%   FP = Fraction of potentiation transitions
%   w  = Weights of states (+/-1)

error(CheckSize(Wp,@ismat));%matrix
error(CheckSize(Wp,@issquare));%square
error(CheckSize(Wm,@(x)samesize(Wp,x),'samesize(Wp)'));%also square matrix of same size
error(CheckSize(fp,@isscalar));
error(CheckValue(fp,@(x) inrange(fp,0,1),'inrange(0,1)'));%fp in [0,1]
error(CheckSize(w,@iscol));
error(CheckValue(w,@(x) all(x.^2==1),'all w = +/-1'));
assert(length(w)==length(Wp));%same size


n=size(Wp,1);
assert(mod(n,2)==0)

q=Wp-Wm;
% w=ones(n,1);
% w(1:(n/2))=-1;

Zinv=ones(size(Wp)) - fp*q - Wm;

area = 2*fp *(1-fp) * sum( (Zinv \ q) * (Zinv \ w) );
% area =  sum( (Zinv \ q) * (Zinv \ w) );


end

