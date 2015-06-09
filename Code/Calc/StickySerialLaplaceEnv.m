function [ As ] = StickySerialLaplaceEnv( s, numstates )
%As=STICKYSERIALLAPLACEENV(s,numstates) maximum Laplace transform of SNR
%curve for sticky serial model
%   As = Laplace transform of SNR curve
%   s         = inverse timescale, Laplace transform parameter
%   numstates = # states

q  = StickySerialLaplaceTrans( s, numstates );
As = StickySerialLaplace( s,numstates,q );


end

