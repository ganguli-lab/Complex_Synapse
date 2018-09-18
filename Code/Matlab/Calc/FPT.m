function [ T ] = FPT( W )
%T=FPT(W) Calculate off diagonal mean first passage times
%   W = transition rates

assert(ismat(W));%matrix
assert(issquare(W));%square

E=ones(size(W));
p=EqProb(W);
Z=inv(E-W);

T=(E*diag(diag(Z)) - Z)/diag(p);

end

