function [ Wp,Wm,w ] = RandTriBuilder( n,sp )
%[Wp,Wm,w]=RANDTRIBUILDER(n,q) build serial model
%   sp = sparsity
%   n  = number of states (even)
%   WP = potentiation transition rates
%   WM = depression transition rates
%   w  = Weights of states (+/-1)

W=RandTrans(n,sp);
w=BinaryWeights(n);

[Wp,Wm,w]=TriangleDcmp(W,0.5,w,NaN);

end

