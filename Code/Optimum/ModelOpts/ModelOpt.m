function [ newWp,newWm,snr,ef ] = ModelOpt( Wp,Wm,tm,varargin)
%[newWp,newWm,snr,ef]=MODELOPT(Wp,Wm,tm) run gradient descent on model to
%maximise snr(t)
%   TM = time value
%   WP = potentiation transition rates
%   WM = depression transition rates
%   snr= snr at T
%   ef = exit flag

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='ModelOpt';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addParameter('UseDerivs',false,@(x) validateattributes(x,{'logical'},{'scalar'},'ModelOpt','UseDerivs'));
    p.addParameter('DispExit',false,@(x) validateattributes(x,{'logical'},{'scalar'},'ModelOpt','DispExit'));
    p.addParameter('TolFun',1e-6,@(x) validateattributes(x,{'numeric'},{'scalar'},'ModelOpt','TolFun'));
    p.addParameter('TolX',1e-10,@(x) validateattributes(x,{'numeric'},{'scalar'},'ModelOpt','TolX'));
    p.addParameter('TolCon',1e-6,@(x) validateattributes(x,{'numeric'},{'scalar'},'ModelOpt','TolCon'));
    p.addParameter('MaxIter',1000,@(x) validateattributes(x,{'numeric'},{'scalar','integer'},'ModelOpt','MaxIter'));
    p.addParameter('Algorithm','interior-point',@(x) parsevalidatestring(x,{'trust-region-reflective','active-set','interior-point','sqp'},'ModelOpt','Algorithm'));
    p.addParameter('Display','off',@(x) parsevalidatestring(x,{'off','iter','iter-detailed','notify','notify-detailed','final','final-detailed'},'ModelOpt','Display'));
    p.addParameter('fp',0.5,@(x) validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1},'ModelOpt','fp'));
end
p.parse(varargin{:});
r=p.Results;

n=length(Wp);
w=BinaryWeights(n);

x0 = Mats2Params(Wp,Wm);            % Starting guess 
lb=zeros(size(x0));
ub=ones(size(x0));
[linconstr_A,linconstr_b]=ParamsConstraints(n);

options = optimset(p.Unmatched,'Algorithm',r.Algorithm,'Display',r.Display,...
    'TolFun', r.TolFun,...  % termination based on function value (of the derivative)
    'TolX', r.TolX,...
    'TolCon',r.TolCon,...
    'MaxIter',r.MaxIter, ...
    'largescale', 'on');

if r.UseDerivs
    options = optimset(options,'GradObj','on');
    [x,snr,ef] = fmincon(@(y)OptFunGrad(y,tm,r.fp,w),x0,...
        linconstr_A,linconstr_b,...
        [],[],lb,ub,[],...
        options);%fun,xo,A,b,Aeq,beq,lb,ub,nonlcon,options
else
    [x,snr,ef] = fmincon(@(y)OptFun(y,tm,r.fp,w),x0,...
        linconstr_A,linconstr_b,...
        [],[],lb,ub,[],...
        options);%fun,xo,A,b,Aeq,beq,lb,ub,nonlcon,options
end
[Wp,Wm]=Params2Mats(x);
snr=-snr;

% [~,~,ix]=SortByEta(0.5*Wp+0.5*Wm,w);
[~,~,ix]=SortByWt(r.fp*Wp+(1-r.fp)*Wm,w,tm);
newWp=Wp(ix,ix);
newWm=Wm(ix,ix);

if r.DispExit
    disp(ExitFlagMsg(ef));
end


end

