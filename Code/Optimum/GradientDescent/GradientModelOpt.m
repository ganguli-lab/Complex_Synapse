function [ newWp,newWm,snr,ef ] = GradientModelOpt( Wp,Wm,tm,varargin)
%[newWp,newWm,snr,ef ]=MODELOPT(Wp,Wm,tm) run gradient descent on model to
%maximise snr(t)
%   TM = time value
%   WP = potentiation transition rates
%   WM = depression transition rates
%   snr= snr at T
%   ef = exit flag

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='GradientModelOpt';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addParameter('DispExit',false,@(x) validateattributes(x,{'logical'},{'scalar'},'GradientModelOpt','DispExit'));
    p.addParameter('TolFun',1e-6,@(x) validateattributes(x,{'numeric'},{'scalar'},'GradientModelOpt','TolFun'));
    p.addParameter('TolX',1e-10,@(x) validateattributes(x,{'numeric'},{'scalar'},'GradientModelOpt','TolFun'));
    p.addParameter('TolCon',1e-6,@(x) validateattributes(x,{'numeric'},{'scalar'},'GradientModelOpt','TolFun'));
    p.addParameter('MaxIter',1000,@(x) validateattributes(x,{'numeric'},{'scalar','integer'},'GradientModelOpt','TolFun'));
    p.addParameter('fp',0.5,@(x) validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1},'GradientModelOpt','fp'));
    p.addParameter('eta',0.5,@(x) validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1},'GradientModelOpt','eta'));
    p.addParameter('w',1,@(x) validateattributes(x,{'numeric'},{'column'},'GradientModelOpt','w'));
end
p.parse(varargin{:});
r=p.Results;

n=length(Wp);

if any(strcmp('w',p.UsingDefaults))
    r.w=BinaryWeights(n);
end


[A,b]=ParamsConstraints(n);

x0 = Mats2Params(Wp,Wm);            % Starting guess 
options = {...
    'TolFun', r.TolFun,...  % termination based on function value (of the derivative)
    'TolX', r.TolX,...
    'TolCon',r.TolCon,...
    'MaxIter',r.MaxIter};

    [x,snr,ef] = GradientDescend(x0,r.eta,A,b,@OptFunGrad,{tm,r.fp,r.w},p.Unmatched,options{:});
[newWp,newWm]=Params2Mats(x);
snr=-snr;
% [~,~,ix]=SortByEta(0.5*Wp+0.5*Wm,w);
% [~,~,ix]=SortByWt(0.5*Wp+0.5*Wm,w,tm);
% newWp=Wp(ix,ix);
% newWm=Wm(ix,ix);


if r.DispExit
    disp(ExitFlagMsg(ef));
end

end

