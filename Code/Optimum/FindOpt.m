function [ Wp,Wm ] = FindOpt( t,n,varargin )
%[Wp,Wm]=FINDOPT(t,n) Find synapse model that maximises SNR(t)
%   t = time value
%   n = #states
%   Wp = potentiation transition rates
%   Wm = depression transition rates

InitRand=true;
UseDerivs=false;
varargin=assignApplicable(varargin);

if InitRand
   [~,~,w]=MakeSMS(ones(1,n-1));
    Wp=RandTrans(n);
    Wm=RandTrans(n);
else
    [Wp,Wm,w]=MakeSMS(ones(1,n-1));
end


[A,b]=ParamsConstraints(n);

x0 = Mats2Params(Wp,Wm);            % Starting guess 
options = optimset('Algorithm','active-set','Display','off');

if UseDerivs
    options = optimset(options,'GradObj','on');
    x = fmincon(@(y)OptFunGrad(y,t,0.5,w),x0,A,b,[],[],[],[],... 
       [],options);
else
    x = fmincon(@(y)OptFun(y,t,0.5,w),x0,A,b,[],[],[],[],... 
       [],options);
end
[Wp,Wm]=Params2Mats(x);

[~,~,ix]=SortByEta(0.5*Wp+0.5*Wm,w);
Wp=Wp(ix,ix);
Wm=Wm(ix,ix);

end

