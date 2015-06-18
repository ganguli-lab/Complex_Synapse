function [ q ] = ShorterSerialLaplaceTrans( s, numstates )
%q=SHORTERSERIALLAPLACETRANS(s,numstates) transition rate for maximum
%Laplace transform of SNR curve for shortened serial model
%   q         = transition rate to end state
%   s         = inverse timescale, Laplace transform parameter
%   numstates = # states

alpha = 1+s + sqrt(s.*(s+2));

q = ((-1)+alpha).^(-1).*((-1).*alpha.^3+alpha.^numstates).^(-1) ...
  .*((-1).*alpha.^3+2.*((-1)+alpha).*alpha.^(1+(1/2).* ...
  numstates)+alpha.^numstates).^(-1).*(alpha.^5.*(1+((-1)+ ...
  alpha).*alpha)+(-2).*alpha.^(4+(1/2).*numstates).*(1+((-1)+ ...
  alpha).*alpha)+((-1)+alpha).*alpha.^(2+numstates).*(1+((-1)+ ...
  alpha).*alpha)+2.*alpha.^(1+(3/2).*numstates).*(1+((-1)+ ...
  alpha).*alpha)+alpha.^(2.*numstates).*((-1)+alpha+(-1).* ...
  alpha.^2)+alpha.^(1+(1/4).*numstates).*((-1).*(1+((-1)+ ...
  alpha).*alpha).*((-1).*alpha+alpha.^numstates).*((-1).* ...
  alpha.^3+alpha.^numstates).*((-2).*alpha.^(1+(1/2).* ...
  numstates).*(2.*alpha+((-1)+alpha).^2.*numstates)+ ...
  alpha.^numstates.*(2.*alpha.^2+numstates+(-1).*alpha.* ...
  numstates)+alpha.^2.*(2+((-1)+alpha).*alpha.*numstates))).^( ...
  1/2));



end
