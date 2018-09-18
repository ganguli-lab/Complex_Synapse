function [ tf ] = isprob( p )
%TF=ISPROB(P) is P a prob measure?
%   is P a row vector of non-neg, summing to 1?
%   TF = logical, [0,1]

tf = isrow(p) && all(p>=0) &&  abs(sum(p)-1) <1e-6;

end

