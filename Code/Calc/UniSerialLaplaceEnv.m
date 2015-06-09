function [ As ] = UniSerialLaplaceEnv( s )
%As=UNISERIALLAPLACEENV(s) maximum Laplace transform of SNR curve for
%uniform serial model
%   As = Laplace transform of SNR curve
%   s         = inverse timescale, Laplace transform parameter

alpha = 1+s + sqrt(s.*(s+2));
xs = 4.50656136546080;


As = ( log(alpha) * (1-xs)^2 ) ./ ( log(xs) * s * (1+xs^2) );


end

