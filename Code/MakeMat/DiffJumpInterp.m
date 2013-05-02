function [ Wp,Wm,w ] = DiffJumpInterp( n,eps )
%[Wp,Wm,w]=DIFFJUMP(n) Possibly better than uniform SMS
%   Detailed explanation goes here

existsAndDefault('eps',1);

q=[ones(1,n/2-1) 0 ones(1,n/2-1)];
q(end)=eps;

[Wp,Wm,w]=MakeSMS(q);

% Wp(n/2,n/2+1)=0;
Wp(n/2,n)=eps;
Wp(n/2,n-1)=1-eps;
% Wm(n/2+1,n/2)=0;
Wm(n/2+1,1)=eps;
Wm(n/2+1,2)=1-eps;

Wp=StochastifyC(Wp);
Wm=StochastifyC(Wm);

end

