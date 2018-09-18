function [ div ] = KLdiv( A,B )
%div=KLDIV(A,B) Kullback-Leibler divergence from A to B
%   A,B = either probability row vectors or Markov matrices

error(CheckSize(A,@(x) isprob(x) || isstochasticD(x),'isprob OR isstochasticD'));
error(CheckSize(B,@(x) isprob(x) || isstochasticD(x),'isprob OR isstochasticD'));
error(CheckSize(B,@(x) samesize(x,A),'samesize(A)'));
% validateattributes(A,{'numeric'},{'matrix','nonnegative'});
% validateattributes(B,{'numeric'},{'size',size(A)});

B=B(A>0);
A=A(A>0);

div=sum(A.*log(A./B));

end

