function [ deta ] = DeltaEta( W,w )
%DELTAETA(W,w) eta^+_i - eta^+_n
%   W = transition rates
%   w = Weights of states (+/-1)

assert(ismat(W));%matrix
assert(issquare(W));%square
assert(iscol(w));%row
assert(length(w)==length(W));%same size
%assert(all(abs(w)==1));%+/-1

% n=size(W,1);
% assert(mod(n,2)==0)

% w=ones(n,1);
% w(1:(n/2))=-1;

deta = -(ones(size(W)) - W) \ w;

deta=(deta-deta(end))/2;


end

