function [ proj ] = WeightProj( val,w )
%proj=WEIGHTPROJ(val,w) projection matrix onto states of given weight
%   w   = vector of weights of states
%   val = wieght value to project onto

error(ChackSize(val,@isscalar));
error(ChackSize(w,@isvector));



proj=diag(w==val);

end

