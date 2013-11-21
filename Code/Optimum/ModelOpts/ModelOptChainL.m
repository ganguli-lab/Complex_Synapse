function [ newqv,A ] = ModelOptChainL( qv,sm,varargin)
%[newqv]=MODELOPTCHAINL(qv,tm) run gradient descent on model
%   sm = Laplace param value
%   qv = nearest neighbour transition rates

UseDerivs=true;
Algorithm='interior-point';
Display='off';
TolFun=1e-6;
TolX=1e-10;
TolCon=1e-6;
MaxIter=1000;
DispExit=false;
varargin=assignApplicable(varargin);

qmin=zeros(size(qv));
qmax=ones(size(qv));

options = optimset('Algorithm',Algorithm,'Display',Display,...
    'TolFun', TolFun,...  % termination based on function value (of the derivative)
    'TolX', TolX,...
    'TolCon',TolCon,...
    'MaxIter',MaxIter, ...
    'largescale', 'on', ...
    varargin{:});

if UseDerivs
    options = optimset(options,'GradObj','on');
    [newqv,A,ef] = fmincon(@(y)OptFunGradChainL(y,sm),qv,...
        [],[],[],[],qmin,qmax,[],... 
        options);
else
    [newqv,A,ef] = fmincon(@(y)OptFunChainL(y,sm),qv,...
        [],[],[],[],qmin,qmax,[],... 
        options);
end
A=-A;

if DispExit
    disp(ExitFlagMsg(ef));
end

end

