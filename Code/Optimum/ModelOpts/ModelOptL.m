function [ newWp,newWm,A,ef ] = ModelOptL( Wp,Wm,sm,varargin)
%[newWp,newWm,A,ef]=MODELOPTL(Wp,Wm,sm) run gradient descent on model to
%maximise A(s)
%   sm = time value
%   WP = potentiation transition rates
%   WM = depression transition rates
%   A  = Laplace Transf value
%   ef = exit flag

UseDerivs=true;
Algorithm='interior-point';
Display='off';
TolFun=1e-6;
TolX=1e-10;
TolCon=1e-6;
MaxIter=1000;
DispExit=false;
fp=0.5;
varargin=assignApplicable(varargin);

n=length(Wp);
w=BinaryWeights(n);

[A,b]=ParamsConstraints(n);

x0 = Mats2Params(Wp,Wm);            % Starting guess 
options = optimset('Algorithm',Algorithm,'Display',Display,...
    'TolFun', TolFun,...  % termination based on function value (of the derivative)
    'TolX', TolX,...
    'TolCon',TolCon,...
    'MaxIter',MaxIter, ...
    'largescale', 'on', ...
    varargin{:});

if UseDerivs
    options = optimset(options,'GradObj','on');
    [x,A,ef] = fmincon(@(y)OptFunGradL(y,sm,fp,w),x0,A,b,...
         [],[],[],[],[],... 
       options);
else
    [x,A,ef] = fmincon(@(y)OptFunL(y,sm,fp,w),x0,A,b,...
         [],[],[],[],[],... 
       options);
end
[Wp,Wm]=Params2Mats(x);

% [~,~,ix]=SortByEta(0.5*Wp+0.5*Wm,w);
[newWp,newWm]=SortByWtEtaS(Wp,Wm,w,fp,sm);
A=-A;

if DispExit
    disp(ExitFlagMsg(ef));
end

end

