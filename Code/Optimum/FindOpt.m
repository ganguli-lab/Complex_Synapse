function [ Wp,Wm ] = FindOpt( t,n )
%FINDOPT Summary of this function goes here
%   Detailed explanation goes here

% [Wp,Wm,w]=MakeSMS(ones(1,n-1));

[~,~,w]=MakeSMS(ones(1,n-1));
Wp=RandTrans(n);
Wm=RandTrans(n);

[A,b]=ParamsConstraints(n);

x0 = Mats2Params(Wp,Wm);            % Starting guess 
options = optimset('Algorithm','active-set','Display','off');

options = optimset(options,'GradObj','on');
x = fmincon(@(y)OptFunGrad(y,t,0.5,w),x0,A,b,[],[],[],[],... 
   [],options);

% x = fmincon(@(y)OptFun(y,t,0.5,w),x0,A,b,[],[],[],[],... 
%    [],options);

[Wp,Wm]=Params2Mats(x);

[~,~,ix]=SortByEta(0.5*Wp+0.5*Wm,w);
Wp=Wp(ix,ix);
Wm=Wm(ix,ix);

end

