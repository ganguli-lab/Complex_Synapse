function [ W ] = AddChain( W,inds,val )
%W=ADDCHAIN(W,INDS,VAL) adds a chain in stochastic matrix, W.
%   chain members link in order
%   INDS = indices of chain members
%   VAL = intra-clique transition rate

assert(ismat(W));
assert(issquare(W));
assert(all(isint(inds)));
assert(isvector(inds));
assert(isscalar(val));


lininds=sub2ind(size(W),inds(1:end-1),inds(2:end));

W(lininds)=val;
W=Stochastify(W);


end

