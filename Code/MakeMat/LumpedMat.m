function [ Wlumped ] = LumpedMat( W,partitions )
%Wtest=LUMPTEST(W,partitions) Trest for lumpability of W wrt partitions
%   Wlumped = lumped version of W
%   W = stochastic matrix /  function (column vec) / measure (row vec)
%   partitions = cell array of vectors containing indices of states in each
%   partition

error(CHeckSize(W,@ismatrix));

[U,V]=LumpProj(partitions);

if iscolumn(W)
    Wlumped = U*W;
elseif isrow(W)
    Wlumped = W*V;
else
    Wlumped = U*W*V;
end

end

