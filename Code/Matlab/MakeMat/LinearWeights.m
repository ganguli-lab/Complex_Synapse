function [ w ] = LinearWeights( n )
%w=LINEARWEIGHTS(n) vector of synaptic weights
%   n=#states

error(CheckSize(n,@isscalar));
error(CheckValue(n,@isint));

w=(-1:2/(n-1):1)';

end

