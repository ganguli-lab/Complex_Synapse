function [ dA ] = DiffArea( Wp,Wm,fp,w,row,col,pm )
%DA=DIFFAREA(WP,WM,FP,w,ROW,COL,PM) derivative of area wrt element
%   WP = potentiation transition rates
%   WM = depression transition rates
%   FP = Fraction of potentiation transitions
%   w = Weights of states (+/-1)
%   ROW,COL = which element we're changing
%   PM = +/-1, if we're changing WP/WM

if isempty(row) || isempty(col)
    dA=0;
    return;
end

assert(ismat(Wp));%matrix
assert(issquare(Wp));%square
assert(samesize(Wp,Wm));%also square matrix of same size
assert(isscalar(fp));
assert(isscalar(row));
assert(isint(row));
assert(isscalar(col));
assert(isint(col));
assert(isscalar(pm));
assert(pm^2==1);

n=size(Wp,1);
assert(mod(n,2)==0)

q=Wp-Wm;
% w=ones(n,1);
% w(1:(n/2))=-1;

pib=ones(1,n);
Zinv=ones(size(Wp)) - fp*q - Wm;

piZ=pib/Zinv;
Zw=Zinv\w;
ZqZ=(Zinv\q)/Zinv;
ZqZw=ZqZ*w;
piZqZ=pib*ZqZ;

fpm=0.5 + pm * (fp-0.5);

dA = fpm * ( piZ(row)*(ZqZw(col)-ZqZw(row)) + piZqZ(row)*(Zw(col)-Zw(row)) ) ...
    + pm * piZ(row)*(Zw(col)-Zw(row));


end

