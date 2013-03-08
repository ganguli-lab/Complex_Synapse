function [ S ] = SNRcurve( t, Wp, Wm, fp, w  )
%S=SNRCURVE(T,WP,WM,FP,w) SNR as function of time
%   T = time values
%   WP = potentiation transition rates
%   WM = depression transition rates
%   FP = Fraction of potentiation transitions
%   w = Weights of states (+/-1)

error(CheckSize(Wp,@ismat));%matrix
error(CheckSize(Wp,@issquare));%square
error(CheckSize(Wm,@(x)samesize(Wp,x),'samesize(Wp)'));%also square matrix of same size
error(CheckSize(fp,@isscalar));
error(CheckValue(fp,@(x) inrange(x,0,1),'inrange(0,1)'));%fp in [0,1]
error(CheckSize(w,@iscol));
error(CheckValue(w,@(x) all(x.^2==1),'all w = +/-1'));


[ qa,ca ] = SpectrumWpm( Wp, Wm, fp, w );

% S = gmdmp(ca.*qa, 1, exp(-outer(qa,t,true)), 1);
S = SNRcurveCaQa(t,ca,qa);

end

