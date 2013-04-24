function [ newWp,newWm ] = ModelOpt( Wp,Wm,tm,varargin)
%[newWp,newWm]=MODELOPT(Wp,Wm,tm) run gradient descent on model
%   T = time value
%   WP = potentiation transition rates
%   WM = depression transition rates

UseDerivs=false;
Algorithm='active-set';
varargin=assignApplicable(varargin);

n=length(Wp);

w=[-ones(n/2,1);ones(n/2,1)];

[A,b]=ParamsConstraints(n);

x0 = Mats2Params(Wp,Wm);            % Starting guess 
options = optimset('Algorithm',Algorithm,'Display','off',varargin{:});

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

