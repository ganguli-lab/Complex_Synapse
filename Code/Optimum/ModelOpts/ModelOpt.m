function [ newWp,newWm,snr,ef ] = ModelOpt( Wp,Wm,tm,varargin)
%[newWp,newWm,snr,ef]=MODELOPT(Wp,Wm,tm) run gradient descent on model to
%maximise snr(t)
%   TM = time value
%   WP = potentiation transition rates
%   WM = depression transition rates
%   snr= snr at T
%   ef = exit flag

UseDerivs=false;
Algorithm='interior-point';
Display='off';
TolFun=1e-6;
TolX=1e-10;
TolCon=1e-6;
MaxIter=1000;
DispExit=false;
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
    [x,snr,ef] = fmincon(@(y)OptFunGrad(y,tm,0.5,w),x0,A,b,...
         [],[],[],[],[],... 
       options);
else
    [x,snr,ef] = fmincon(@(y)OptFun(y,tm,0.5,w),x0,A,b,...
         [],[],[],[],[],... 
       options);
end
[Wp,Wm]=Params2Mats(x);
snr=-snr;

% [~,~,ix]=SortByEta(0.5*Wp+0.5*Wm,w);
[~,~,ix]=SortByWt(0.5*Wp+0.5*Wm,w,tm);
newWp=Wp(ix,ix);
newWm=Wm(ix,ix);

if DispExit
    disp(ExitFlagMsg(ef));
end


end

