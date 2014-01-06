function [ Wp,Wm,w ] = NonuniBuilder( n,x )
%[Wp,Wm,w]=NonuniBuilder(n,x) build non-uniform multistate model
%   x  = ratio of nearby transition rates
%   n  = number of states (even)
%   WP = potentiation transition rates
%   WM = depression transition rates
%   w  = Weights of states (+/-1)


[Wp,Wm]=CascadeOriginal(x,x,n,0);
Wp=(1-x)*Wp;
Wm=(1-x)*Wm;
w=LinearWeights(n);

end

