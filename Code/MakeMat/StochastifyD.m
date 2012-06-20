function [ newW ] = StochastifyD( W )
%NEWW=STOCHASTIFY(W) turn into discrete time stochastic matrix
%   Scales rows so that sum_j W_ij = 1

assert(ismat(W));
assert(issquare(W));

newW=diag(1./sum(W,2))*W;


end

