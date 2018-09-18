function [ c ] = AreaCoeffSorted( Wp, Wm, fp,w )
%C=AREACOEFF(WP,WM,FP,w) coeff of p^\infty_k w_k in area formula
%   WP = potentiation transition rates
%   WM = depression transition rates
%   FP = Fraction of potentiation transitions
%   w  = Weights of states (+/-1)

assert(ismat(Wp));%matrix
assert(issquare(Wp));%square
assert(samesize(Wp,Wm));%also square matrix of same size
assert(isscalar(fp));
assert(0<=fp && fp<=1);%fp in [0,1]
assert(iscol(w));%row
assert(length(w)==length(Wp));%same size
assert(all(abs(w)==1));%+/-1


c=AreaCoeff(Wp,Wm,fp);
deta=DeltaEta(fp*Wp+(1-fp)*Wm,w);

[~,ix]=sort(deta,'descend');
c=c(ix);

end

