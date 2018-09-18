function [ W ] = RandTrans( len, sparsity )
%W=RANDTRANS(LEN) Random Markov matrix (cts time)
%   LEN = # states

error(CheckSize(len,@isscalar));
error(CheckValue(len,@isint));


W=StochastifyC(StochastifyD(rand(len)));

if exist('sparsity','var')
    S=rand(size(W));
    S=(S<=sparsity);
    W=StochastifyC(W.*S);
end%if sparsity exists

end

