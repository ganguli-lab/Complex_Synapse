function [ As ] = ShorterSerialLaplaceEnv( s, numstates )
%As=SHORTERSERIALLAPLACEENV(s,numstates) maximum Laplace transform of SNR
%curve for shortened serial model
%   As = Laplace transform of SNR curve
%   s         = inverse timescale, Laplace transform parameter
%   numstates = # states

q  = ShorterSerialLaplaceTrans( s, numstates );
As = ShorterSerialLaplace( s,numstates,q );


end

