function [ c ] = AreaCoeff( Wp, Wm, fp )
%C=AREACOEFF(WP,WM,FP) coeff of p^\infty_k w_k in area formula
%   WP = potentiation transition rates
%   WM = depression transition rates
%   FP = Fraction of potentiation transitions

assert(ismat(Wp));%matrix
assert(issquare(Wp));%square
assert(samesize(Wp,Wm));%also square matrix of same size
assert(isscalar(fp));
assert(0<=fp && fp<=1);%fp in [0,1]

q=Wp-Wm;
pib=ones(1,length(Wp));
Zinv=ones(size(Wp)) - fp*q - Wm;

piZ=pib/Zinv;
piZqZ=(piZ*q)/Zinv;

c=piZqZ./piZ;


end

