function [ Aenv,sbnd] = SerialLaplaceEnv( s,numstates )
%As=SerialLaplaceEnv(s,numstates) maximum Laplace transform of SNR
%curve for serial topology models
%   s         = inverse timescale, Laplace transform parameter
%   numstates = # states
%   Aenv = max Laplace transform of SNR curve
%   sbnd = boundaries of regimes [sticky-uniform-shorter-binary]

Aenv=UniSerialLaplaceEnv(s);

AenvUni=UniSerialLaplace(s,numstates);
AenvBinary=UniSerialLaplace(s,2);
AenvSticky=StickySerialLaplaceEnv(s,12);

stickymax=StickySerialLaplaceEnvMax(numstates);
[smin,smax]=UniSerialLaplaceEnvMin(numstates);


Aenv(s>smax)=AenvBinary(s>smax);
Aenv(s<smin)=AenvUni(s<smin);
Aenv(s<stickymax)=AenvSticky(s<stickymax);


sbnd=[stickymax smin smax];

end

