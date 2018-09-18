function [ W ] = AddClique( W,inds,val )
%W=ADDCLIQUE(W,INDS,VAL) adds an all-to-all clique in stochastic matrix, W.
%   INDS = indices of clique members
%   VAL = intra-clique transition rate

assert(ismat(W));
assert(issquare(W));
assert(all(isint(inds)));
assert(isvector(inds));
assert(isscalar(val));


% W(inds,inds)=val;
W(inds,inds)=val*rand(length(inds));
W=Stochastify(W);


end

