function [ newWp,newWm,snr,ef ] = GradientModelOpt( Wp,Wm,tm,varargin)
%[newWp,newWm,snr,ef ]=MODELOPT(Wp,Wm,tm) run gradient descent on model to
%maximise snr(t)
%   TM = time value
%   WP = potentiation transition rates
%   WM = depression transition rates
%   snr= snr at T
%   ef = exit flag

eta=1e-3;
TolFun=1e-6;
TolX=1e-10;
TolCon=1e-6;
MaxIter=1000;
DispExit=false;
n=length(Wp);
w=BinaryWeights(n);
fp=0.5;
varargin=assignApplicable(varargin);


[A,b]=ParamsConstraints(n);

x0 = Mats2Params(Wp,Wm);            % Starting guess 
options = {...
    'TolFun', TolFun,...  % termination based on function value (of the derivative)
    'TolX', TolX,...
    'TolCon',TolCon,...
    'MaxIter',MaxIter, ...
    varargin{:}};

    [x,snr,ef] = GradientDescend(x0,eta,A,b,@OptFunGrad,tm,fp,w,options{:});
[newWp,newWm]=Params2Mats(x);
snr=-snr;
% [~,~,ix]=SortByEta(0.5*Wp+0.5*Wm,w);
% [~,~,ix]=SortByWt(0.5*Wp+0.5*Wm,w,tm);
% newWp=Wp(ix,ix);
% newWm=Wm(ix,ix);


if DispExit
    disp(ExitFlagMsg(ef));
end

end

