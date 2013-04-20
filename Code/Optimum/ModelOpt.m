function [ newWp,newWm ] = ModelOpt( Wp,Wm,tm,varargin)
%MODELOPT Summary of this function goes here
%   Detailed explanation goes here

UseDerivs=false;
varargin=assignApplicable(varargin);

n=length(Wp);

w=[-ones(n/2,1);ones(n/2,1)];

[A,b]=ParamsConstraints(n);

x0 = Mats2Params(Wp,Wm);            % Starting guess 
options = optimset('Algorithm','active-set','Display','off');

if UseDerivs
    options = optimset(options,'GradObj','on');
    x = fmincon(@(y)OptFunGrad(y,tm,0.5,w),x0,A,b,[],[],[],[],... 
       [],options);
else
    x = fmincon(@(y)OptFun(y,tm,0.5,w),x0,A,b,[],[],[],[],... 
       [],options);
end
[Wp,Wm]=Params2Mats(x);

[~,~,ix]=SortByEta(0.5*Wp+0.5*Wm,w);
newWp=Wp(ix,ix);
newWm=Wm(ix,ix);

end

