function [ newW ] = Stochastify( W )
%NEWW=STOCHASTIFY(W) turn into discrete time cts time stochastic matrix
%with cts time 
%   Changes diagonal terms to W_ii = sum_{j~=i} W_ij

assert(ismat(W));
assert(issquare(W));

newW=StochastifyC(StochastifyD(W));


end

