function [ As ] = UniSerialLaplace( s, numstates )
%s=UNISERIALLAPLACE(s,numstates) Laplace transform of SNR curve for uniform
%serial model
%   As = Laplace transform of SNR curve
%   s         = inverse timescale, Laplace transform parameter
%   numstates = # states


alpha = 1+s + sqrt(s.*(s+2));
x = alpha.^(numstates/2);


As = ( 2 * (1-x).^2 ) ./ ( numstates * s .* (1+x.^2) );





end

