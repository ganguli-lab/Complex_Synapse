function [ newqv,A,ef ] = ModelOptChainHomL( qv,sm,varargin)
%[newqv,A,ef]=MODELOPTCHAINHOML(qv,sm) run gradient descent on model to
%maximise A(s)
%   sm = Laplace param value
%   qv = nearest neighbour transition rates
%   A  = Laplace Transf value
%   ef = exit flag

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='ModelOptChainHomL';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addParameter('UseDerivs',true,@(x) validateattributes(x,{'logical'},{'scalar'},'ModelOptChainHomL','UseDerivs'));
    p.addParameter('DispExit',false,@(x) validateattributes(x,{'logical'},{'scalar'},'ModelOptChainHomL','DispExit'));
    p.addParameter('TolFun',1e-6,@(x) validateattributes(x,{'numeric'},{'scalar'},'ModelOptChainHomL','TolFun'));
    p.addParameter('TolX',1e-10,@(x) validateattributes(x,{'numeric'},{'scalar'},'ModelOptChainHomL','TolFun'));
    p.addParameter('TolCon',1e-6,@(x) validateattributes(x,{'numeric'},{'scalar'},'ModelOptChainHomL','TolFun'));
    p.addParameter('MaxIter',1000,@(x) validateattributes(x,{'numeric'},{'scalar','integer'},'ModelOptChainHomL','TolFun'));
    p.addParameter('Algorithm','interior-point',@(x) validatestring(x,{'trust-region-reflective','active-set','interior-point','sqp'},'ModelOptChainHomL','TolFun'));
    p.addParameter('Display','off',@(x) validatestring(x,{'off','iter','iter-detailed','notify','notify-detailed','final','final-detailed'},'ModelOptChainHomL','TolFun'));
    p.addParameter('fp',0.5,@(x) validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1},'ModelOptChainHomL','fp'));
end
p.parse(varargin{:});
r=p.Results;

qmin=zeros(size(qv));
qmax=ones(size(qv));
qmax(length(qmax)/2+1:end)=Inf;

options = optimset(p.Unmatched,'Algorithm',r.Algorithm,'Display',r.Display,...
    'TolFun', r.TolFun,...  % termination based on function value (of the derivative)
    'TolX', r.TolX,...
    'TolCon',r.TolCon,...
    'MaxIter',r.MaxIter, ...
    'largescale', 'on');

if r.UseDerivs
    options = optimset(options,'GradObj','on');
    [newqv,A,ef] = fmincon(@(y)OptFunGradChainHomL(y,sm,r.fp),qv,...
        [],[],[],[],qmin,qmax,[],... 
        options);
else
    [newqv,A,ef] = fmincon(@(y)OptFunChainHomL(y,sm,r.fp),qv,...
        [],[],[],[],qmin,qmax,[],... 
        options);
end
A=-A;

if r.DispExit
    disp(ExitFlagMsg(ef));
end

end

