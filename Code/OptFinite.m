function [ dSdWp,dSdWm ] = OptFinite( n )
%OPTFINITE Summary of this function goes here
%   Detailed explanation goes here
eps=0.001;

tm=((n-1)/3)^2;

qq=ones(1,n-1);
% qq=ones(1,n+1);
% qq(end)=eps;
[Wp,Wm,w]=MakeSMS(qq);

[~,dSdWp,dSdWm]=GradSNR(tm,Wp,Wm,0.5,w);

end

