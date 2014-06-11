function [ tf ] = isstochasticD( W )
%tf=ISSTOCHASTIC(W) is it a valid discrete time stochastic matrix, i.e. do rows
%sum to one, and are all elements positive?
%   tf=logical[0,1]

tf = ismat(W) && issquare(W) && all(abs(sum(W,2)-1)<1e-7) && all(all(W>=-1e-7));


end

