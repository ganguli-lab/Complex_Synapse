function [ tf ] = isstochastic( W )
%tf=ISSTOCHASTIC(W) is it a valid continuous time stochastic matrix, i.e.
%do rows sum to zero, and are off diagonal elements positive?
%   tf=logical[0,1]

tf = ismat(W) && issquare(W) && all(abs(sum(W,2))<1e-7) && all(all(W-diag(diag(W)) >=0));


end

