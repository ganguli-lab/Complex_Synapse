function [ As ] = StickySerialLaplace( s,numstates,q )
%As=STICKYSERIALLAPLACE(s,numstates,q) Laplace transform of SNR curve for
%sticky serial model
%   As = Laplace transform of SNR curve
%   s         = inverse timescale, Laplace transform parameter
%   numstates = # states
%   q         = transition rate from end state

alpha = 1+s + sqrt(s.*(s+2));


As=2*q./((2+(numstates-2)*q).*s)...
    .*( 1 - 2*q.*alpha.^(numstates/2)...
      ./( 1 + alpha.^numstates + (q-1).*alpha .*(1+alpha.^(numstates-2)) ) );

end

