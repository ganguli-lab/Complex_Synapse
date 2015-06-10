function [ As ] = ShorterSerialLaplace( s,numstates,q )
%As=SHORTERSERIALLAPLACE(s,numstates,q) Laplace transform of SNR curve for
%shortened serial model
%   As = Laplace transform of SNR curve
%   s         = inverse timescale, Laplace transform parameter
%   numstates = # states
%   q         = transition rate to end state

alpha = 1+s + sqrt(s.*(s+2));


As=2.*((-2)+numstates+2.*q).^(-1).*(1+2.*alpha.^(1+(1/2).* ...
  numstates).*((-1)+alpha+(-1).*alpha.^2+((-1)+alpha).^2.*q).* ...
  (alpha.^2+alpha.^3.*((-1)+q)+(-1).*alpha.^4.*((-1)+q)+ ...
  alpha.^numstates.*(1+((-1)+alpha).*alpha+((-1)+alpha).*q)) ...
  .^(-1)).*s.^(-1);

end

