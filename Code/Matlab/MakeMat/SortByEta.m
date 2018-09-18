function [ newW,neww,ix ] = SortByEta( W,w )
%[newW,neww,ix]=SORTBYETA(W,w) Put states in order of decreasing eta^+
%   W = transition rates
%   w = Weights of states (+/-1)
%   ix=sort order
%   newW=W(ix,ix)
%   neww=w(ix)

assert(ismat(W));%matrix
assert(issquare(W));%square
assert(iscol(w));%row
assert(length(w)==length(W));%same size
assert(all(abs(w)==1));%+/-1

deta=DeltaEta(W,w);

[~,ix]=sort(deta,'descend');

newW=W(ix,ix);
neww=w(ix);


end

