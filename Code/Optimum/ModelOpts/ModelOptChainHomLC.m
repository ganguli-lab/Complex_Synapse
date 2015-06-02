function [ newqv,A,ef ] = ModelOptChainHomLC( qv,sm,varargin)
%[newqv,A,ef]=MODELOPTCHAINHOMLC(qv,sm) run gradient descent on model to
%maximise A(s)
%   sm = Laplace param value
%   qv = nearest neighbour transition rates
%   A  = Laplace Transf value
%   ef = exit flag

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='ModelOptChainHomLC';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addParameter('UseDerivs',true,@(x) validateattributes(x,{'logical'},{'scalar'},'ModelOptChainHomLC','UseDerivs'));
    p.addParameter('DispExit',false,@(x) validateattributes(x,{'logical'},{'scalar'},'ModelOptChainHomLC','DispExit'));
    p.addParameter('TolFun',1e-6,@(x) validateattributes(x,{'numeric'},{'scalar'},'ModelOptChainHomLC','TolFun'));
    p.addParameter('TolX',1e-10,@(x) validateattributes(x,{'numeric'},{'scalar'},'ModelOptChainHomLC','TolFun'));
    p.addParameter('TolCon',1e-6,@(x) validateattributes(x,{'numeric'},{'scalar'},'ModelOptChainHomLC','TolFun'));
    p.addParameter('MaxIter',1000,@(x) validateattributes(x,{'numeric'},{'scalar','integer'},'ModelOptChainHomLC','TolFun'));
    p.addParameter('Algorithm','interior-point',@(x) validatestring(x,{'trust-region-reflective','active-set','interior-point','sqp'},'ModelOptChainHomLC','TolFun'));
    p.addParameter('Display','off',@(x) validatestring(x,{'off','iter','iter-detailed','notify','notify-detailed','final','final-detailed'},'ModelOptChainHomLC','TolFun'));
    p.addParameter('fp',0.5,@(x) validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1},'ModelOptChainHomLC','fp'));
end
p.parse(varargin{:});
r=p.Results;

[ Ac,b ] = HomChainConstr( length(qv)/4+1, r.fp );


options = optimset(p.Unmatched,'Algorithm',r.Algorithm,'Display',r.Display,...
    'TolFun', r.TolFun,...  % termination based on function value (of the derivative)
    'TolX', r.TolX,...
    'TolCon',r.TolCon,...
    'MaxIter',r.MaxIter, ...
    'largescale', 'on');

if r.UseDerivs
    options = optimset(options,'GradObj','on');
    [newqv,A,ef] = fmincon(@(y)OptFunGradChainHomLC(y,sm,r.fp),qv,...
        Ac,b,[],[],[],[],[],... 
        options);
else
    [newqv,A,ef] = fmincon(@(y)OptFunChainHomLC(y,sm,r.fp),qv,...
        Ac,b,[],[],[],[],[],... 
        options);
end
A=-A;

if any(Ac*newqv'>b)
    A=NaN;
end


if r.DispExit
    disp(ExitFlagMsg(ef));
end

end

