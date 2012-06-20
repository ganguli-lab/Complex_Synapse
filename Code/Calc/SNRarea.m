function [ area ] = SNRarea( Wp, Wm, fp, w )
%AREA=SNRAREA(WP,WM,FP,w) Area under SNR curve for complex synapse (cts time)
%   WP = potentiation transition rates
%   WM = depression transition rates
%   FP = Fraction of potentiation transitions
%   w  = Weights of states (+/-1)

assert(ismat(Wp));%matrix
assert(issquare(Wp));%square
assert(samesize(Wp,Wm));%also square matrix of same size
assert(isscalar(fp));
assert(0<=fp && fp<=1);%fp in [0,1]
assert(iscol(w));%row
assert(length(w)==length(Wp));%same size
assert(all(abs(w)==1));%+/-1


n=size(Wp,1);
assert(mod(n,2)==0)

q=Wp-Wm;
% w=ones(n,1);
% w(1:(n/2))=-1;

Zinv=ones(size(Wp)) - fp*q - Wm;

area = 2*fp *(1-fp) * sum( (Zinv \ q) * (Zinv \ w) );
% area =  sum( (Zinv \ q) * (Zinv \ w) );


end

