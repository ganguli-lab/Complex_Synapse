function [ w ] = BinaryWeights( n )
%w=BINARYWEIGHTS(n) vector of synaptic weights
%   n=#states

error(CheckSize(n,@isscalar));
error(CheckValue(n,@isint));

w=[-ones(n/2,1);ones(n/2,1)];

end

