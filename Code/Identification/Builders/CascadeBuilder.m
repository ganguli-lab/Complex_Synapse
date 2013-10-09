function [ Wp,Wm,w ] = CascadeBuilder( n,x )
%[Wp,Wm,w]=CascadeBuilder(n,x) build cascade model
%   x  = ratio of nearby transition rates
%   n  = number of states (even)
%   WP = potentiation transition rates
%   WM = depression transition rates
%   w  = Weights of states (+/-1)


[Wp,Wm,w]=CascadeOriginal(x,x,n);

end

