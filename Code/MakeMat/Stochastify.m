function [ newW ] = Stochastify( W )
%NEWW=STOCHASTIFY(W) turn into cts time stochastic matrix
%   Changes diagonal terms to W_ii = sum_{j~=i} W_ij

assert(ismat(W));
assert(issquare(W));

newW=W-diag(sum(W,2));


end

