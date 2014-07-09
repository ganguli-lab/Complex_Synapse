function [ newWp,newWm,A,ef ] = ModelOptL( Wp,Wm,sm,varargin)
%[newWp,newWm,A,ef]=MODELOPTL(Wp,Wm,sm) run gradient descent on model to
%maximise A(s)
%   sm = time value
%   WP = potentiation transition rates
%   WM = depression transition rates
%   A  = Laplace Transf value
%   ef = exit flag

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='ModelOptL';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addParameter('UseDerivs',true,@(x) validateattributes(x,{'logical'},{'scalar'},'ModelOptL','UseDerivs'));
    p.addParameter('DispExit',false,@(x) validateattributes(x,{'logical'},{'scalar'},'ModelOptL','DispExit'));
    p.addParameter('TolFun',1e-6,@(x) validateattributes(x,{'numeric'},{'scalar'},'ModelOptL','TolFun'));
    p.addParameter('TolX',1e-10,@(x) validateattributes(x,{'numeric'},{'scalar'},'ModelOptL','TolFun'));
    p.addParameter('TolCon',1e-6,@(x) validateattributes(x,{'numeric'},{'scalar'},'ModelOptL','TolFun'));
    p.addParameter('MaxIter',1000,@(x) validateattributes(x,{'numeric'},{'scalar','integer'},'ModelOptL','TolFun'));
    p.addParameter('Algorithm','interior-point',@(x) validatestring(x,{'trust-region-reflective','active-set','interior-point','sqp'},'ModelOptL','TolFun'));
    p.addParameter('Display','off',@(x) validatestring(x,{'off','iter','iter-detailed','notify','notify-detailed','final','final-detailed'},'ModelOptL','TolFun'));
    p.addParameter('fp',0.5,@(x) validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1},'ModelOptL','fp'));
end
p.parse(varargin{:});
r=p.Results;

fp=r.fp;

n=length(Wp);
w=BinaryWeights(n);

[A,b]=ParamsConstraints(n);

x0 = Mats2Params(Wp,Wm);            % Starting guess 
options = optimset(p.Unmatched,'Algorithm',r.Algorithm,'Display',r.Display,...
    'TolFun', r.TolFun,...  % termination based on function value (of the derivative)
    'TolX', r.TolX,...
    'TolCon',r.TolCon,...
    'MaxIter',r.MaxIter, ...
    'largescale', 'on');

if r.UseDerivs
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

if r.DispExit
    disp(ExitFlagMsg(ef));
end

end

