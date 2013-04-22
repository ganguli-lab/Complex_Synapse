function [ S ] = SNRcurveWq( t, W, q ,fp, w, varargin  )
%S=SNRCURVE(T,W,q,w) SNR as function of time
%   T = time values
%   W = forgetting transition rates: f^+W^+ + f^-W^-
%   q = encoding transition rates: W^+ - W^-
%   w = Weights of states (+/-1)

UseExpM=false;
varargin=assignApplicable(varargin);

if UseExpM
    S=zeros(size(t));
    p=EqProb(W);
    for i=1:numel(t)
        S(i) = p*q*expm(W*t(i))*w;
    end
    S=2*fp*(1-fp)*S;
else
    [ qa,ca ] = SpectrumWq( W,q ,fp, w );
    S = SNRcurveCaQa(t,ca,qa);
end

end

