function [ Wp,Wm,w ] = DiffJump( n,eps,q1,q2,q3 )
%[Wp,Wm,w]=DIFFJUMP(n) Possibly better than uniform SMS
%   Detailed explanation goes here

existsAndDefault('eps',1);
existsAndDefault('q1',1);
existsAndDefault('q2',1);
existsAndDefault('q3',1);

q=[q1*ones(1,n/2-1) 0 q2*ones(1,n/2-1)];
q(1)=eps;

[Wp,Wm,w]=MakeSMS(q);

% Wp(n/2,n/2+1)=0;
Wp(n/2,n)=q3;
% Wm(n/2+1,n/2)=0;
Wm(n/2+1,1)=q3;

Wp=StochastifyC(Wp);
Wm=StochastifyC(Wm);

end

