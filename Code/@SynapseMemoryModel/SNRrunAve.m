function [ A ] = SNRrunAve( obj,t)
%A=SynapseMemoryModel.SNRRUNAVE(tau) running maverage of SNR
%curve for complex synapse (cts time), i.e. (Laplace Transform @ s=1/t)/t.
%   A(tau) = int exp(-t/tau)*SNR(t) dt / tau
%   tau  = time = 1/parameter of Laplace transform
%   WP = potentiation transition rates
%   WM = depression transition rates
%   FP = Fraction of potentiation transitions
%   w  = Weights of states (+/-1)
A = obj.SNRlaplace(1./t)./t;


end

