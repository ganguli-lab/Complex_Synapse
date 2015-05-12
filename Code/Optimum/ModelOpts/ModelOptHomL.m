function [ newWp,newWm,newQ,A,ef ] = ModelOptHomL( Wp,Wm,Q,sm,varargin)
%[newWp,newWm,newQ,A,ef]=MODELOPTHomL(Wp,Wm,Q,sm) run gradient descent on model to
%maximise A(s)
%   sm = time value
%   WP = potentiation transition rates
%   WM = depression transition rates
%   Q  = activity independent transition rates
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

[Ac,b]=ParamsConstraintsHom(n);

x0 = Mats2ParamsHom(Wp,Wm,Q);            % Starting guess 
options = optimset(p.Unmatched,'Algorithm',r.Algorithm,'Display',r.Display,...
    'TolFun', r.TolFun,...  % termination based on function value (of the derivative)
    'TolX', r.TolX,...
    'TolCon',r.TolCon,...
    'MaxIter',r.MaxIter, ...
    'largescale', 'on');

if r.UseDerivs
    options = optimset(options,'GradObj','on');
    [x,A,ef] = fmincon(@(y)OptFunGradHomL(y,sm,fp,w),x0,Ac,b,...
         [],[],[],[],[],... 
       options);
else
    [x,A,ef] = fmincon(@(y)OptFunHomL(y,sm,fp,w),x0,Ac,b,...
         [],[],[],[],[],... 
       options);
end
[Wp,Wm,Q]=Params2MatsHom(x);

% [~,~,ix]=SortByEta(0.5*Wp+0.5*Wm,w);
[newWp,newWm,newQ]=SortByWtEtaSHom(Wp,Wm,Q,w,fp,sm);
A=-A;

if any(Ac*x>b)
    A=NaN;
end


if r.DispExit
    disp(ExitFlagMsg(ef));
end

end

