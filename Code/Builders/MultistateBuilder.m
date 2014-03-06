function [ Wp,Wm,w ] = MultistateBuilder( n,q )
%[Wp,Wm,w]=SERIALBUILDER(n,q) build serial model
%   q  = transition rates
%   n  = number of states (even)
%   WP = potentiation transition rates
%   WM = depression transition rates
%   w  = Weights of states (+/-1)


[Wp,Wm]=MakeSMS(q*ones(1,n-1));
w=LinearWeights(n);
end

