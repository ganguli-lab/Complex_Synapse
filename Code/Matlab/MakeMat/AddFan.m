function [ W ] = AddFan( W,indsfrom,indsto,val )
%W=ADDFAN(W,INDSFROM,INDSTO,VAL) adds a fan in stochastic matrix, W.
%   every member of start links to every element of end
%   INDSFROM = indices of start of fan
%   INDSTO = indices of end of fan
%   VAL = intra-clique transition rate

assert(ismat(W));
assert(issquare(W));
assert(all(isint(indsfrom)));
assert(isvector(indsfrom));
assert(all(isint(indsto)));
assert(isvector(indsto));
assert(isscalar(val));


W(indsfrom,indsto)=val;
W=Stochastify(W);


end

