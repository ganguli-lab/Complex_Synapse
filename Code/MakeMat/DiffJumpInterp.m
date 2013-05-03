function [ Wp,Wm,w ] = DiffJumpInterp( n,eps )
%[Wp,Wm,w]=DIFFJUMP(n) Possibly better than uniform SMS
%   Detailed explanation goes here

existsAndDefault('eps',1);

q=[ones(1,n/2-1) 0 ones(1,n/2-1)];
%
% q(end)=eps
%
% q(n/2-1)=eps;
%
q(n/2+1)=eps;
q(1)=eps;

[Wp,~,w]=MakeSMS(q);

%Wp(n/2,n-1)=1-eps
%
% Wp(n/2,n)=eps;
%
% Wp(n/2-1,n/2)=eps;
% Wp(n/2-1,n)=1-eps;
%
Wp(n/2,n)=1;
Wp(n/2+1,n/2+3)=1-eps;
Wp(1,3)=1-eps;

Wp=StochastifyC(Wp);

Wm=rot90(Wp,2);

end

