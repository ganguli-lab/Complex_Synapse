function [ val ] = crq( Wp,Wm,fp,w )
%val=CRQ(WP,WM,FP,w) sum_a I_a sqrt(tau_a)
%   WP = potentiation transition rates
%   WM = depression transition rates
%   FP = Fraction of potentiation transitions
%   w = Weights of states (+/-1)

[qa,ca]=SpectrumWpm(Wp,Wm,fp,w);
val=ca'*sqrt(qa);

end

