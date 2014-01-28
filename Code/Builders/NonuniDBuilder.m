function [ Wp,Wm,w ] = NonuniDBuilder( n,x )
%[Wp,Wm,w]=NonuniBuilder(n,x) build non-uniform multistate model
%   x  = ratio of nearby transition rates
%   n  = number of states (even)
%   WP = potentiation transition rates
%   WM = depression transition rates
%   w  = Weights of states (+/-1)

% qp=(x.^abs((1:n-1)-n/2));
qp=x*ones(1,n-1);
qm=x.^abs((1:n-1)-n/2);

[Wp,Wm]=MakeMultistate(qp,qm);
w=LinearWeights(n);

end

