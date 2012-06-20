function [ S ] = SNRcurveWq( t, W, q ,fp, w  )
%S=SNRCURVE(T,W,q,w) SNR as function of time
%   T = time values
%   W = forgetting transition rates: f^+W^+ + f^-W^-
%   q = encoding transition rates: W^+ - W^-
%   w = Weights of states (+/-1)



[ qa,ca ] = SpectrumWq( W,q ,fp, w );

% S = gmdmp(ca.*qa, 1, exp(-outer(qa,t,true)), 1);
S = SNRcurveCaQa(t,ca,qa);

end

