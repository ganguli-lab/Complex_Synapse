function [ Wp,Wm,w ] = DiffJump( n )
%[Wp,Wm,w]=DIFFJUMP(n) Possibly better than uniform SMS
%   Detailed explanation goes here
[Wp,Wm,w]=MakeSMS(ones(1,n-1));

Wp(n/2,n/2+1)=0;
Wp(n/2,n)=1;
Wm(n/2+1,n/2)=0;
Wm(n/2+1,1)=1;

end

