function [ S ] = SNRcurve( t, Wp, Wm, fp, w  )
%S=SNRCURVE(T,WP,WM,FP,w) SNR as function of time
%   T = time values
%   WP = potentiation transition rates
%   WM = depression transition rates
%   FP = Fraction of potentiation transitions
%   w = Weights of states (+/-1)

[ qa,ca ] = Spectrum( Wp, Wm, fp, w );

% S = gmdmp(ca.*qa, 1, exp(-outer(qa,t,true)), 1);
S = SNRcurveCaQa(t,ca,qa);

end

