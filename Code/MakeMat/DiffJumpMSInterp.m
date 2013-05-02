function [ Wp,Wm,w ] = DiffJumpMSInterp( n,eps )
%[Wp,Wm,w]=DIFFJUMP(n) Possibly better than uniform SMS
%   Detailed explanation goes here

existsAndDefault('eps',1);

q=[ones(1,n/2-1) 1-eps ones(1,n/2-1)];

[Wp,Wm,w]=MakeSMS(q);

% Wp(n/2,n/2+1)=0;
Wp(n/2,n)=eps;
% Wm(n/2+1,n/2)=0;
Wm(n/2+1,1)=eps;

Wp=StochastifyC(Wp);
Wm=StochastifyC(Wm);

end

