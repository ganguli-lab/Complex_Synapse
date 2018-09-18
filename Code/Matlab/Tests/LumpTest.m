function [ Wtest ] = LumpTest( W,partitions )
%Wtest=LUMPTEST(W,partitions) Test for lumpability of W wrt partitions
%   Wtest = quantity that should be all zeros
%   W = stochastic matrix /  function (column vector)
%   partitions = cell array of vectors containing indices of states in each
%                partition

error(CheckSize(W,@ismatrix));

[U,V]=LumpProj(partitions);

if iscol(W)
    Wtest = V*U*W-W;
else
    Wtest = V*U*W*V-W*V;
end

end

