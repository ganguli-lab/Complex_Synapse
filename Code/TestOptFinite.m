function [ KTmultsp,KTmultsm,fliptest ] = TestOptFinite( n )
%OPTFINITE Summary of this function goes here
%   Detailed explanation goes here

[dSdWp,dSdWm]=OptFinite(n);

KTdiag=diag(dSdWp,1);
KTmultsp=[KTdiag;0]*ones(1,n)-dSdWp;

KTdiag=diag(dSdWm,-1);
KTmultsm=[0;KTdiag]*ones(1,n)-dSdWm;

fliptest=dSdWp-rot90(dSdWm,2);

end

