function [ Wp,Wm,w ] = SerialJumpBuilder( n,q )
%[Wp,Wm,w]=SERIALJUMPBUILDER(n,q) build serial-jump model
%   q  = transition rates
%   n  = number of states (even)
%   WP = potentiation transition rates
%   WM = depression transition rates
%   w  = Weights of states (+/-1)


[Wp,Wm,w]=DiffJump(n,q,q,q,q);

end

