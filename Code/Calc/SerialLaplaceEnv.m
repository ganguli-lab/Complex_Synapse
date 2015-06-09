function [ Aenv ] = SerialLaplaceEnv( s,numstates )
%As=SerialLaplaceEnv(s,numstates) maximum Laplace transform of SNR
%curve for serial topology models
%   Aenv = max Laplace transform of SNR curve
%   s         = inverse timescale, Laplace transform parameter
%   numstates = # states

Aenv=UniSerialLaplaceEnv(s);

AenvUni=UniSerialLaplace(s,numstates);
AenvBinary=UniSerialLaplace(s,2);
AenvSticky=StickySerialLaplaceEnv(s,12);

stickymax=StickySerialLaplaceEnvMax(numstates);
[smin,smax]=UniSerialLaplaceEnvMin(numstates);


Aenv(s>smax)=AenvBinary(s>smax);
Aenv(s<smin)=AenvUni(s<smin);
Aenv(s<stickymax)=AenvSticky(s<stickymax);

end

