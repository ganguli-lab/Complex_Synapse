function [ newWp,newWm,A,ef ] = ModelOptL( Wp,Wm,sm,varargin)
%[newWp,newWm,A,ef]=MODELOPTL(Wp,Wm,sm) run gradient descent on model to
%maximise A(sm)
%   sm = inverse time value
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
    p.addParameter('TolX',1e-10,@(x) validateattributes(x,{'numeric'},{'scalar'},'ModelOptL','TolX'));
    p.addParameter('TolCon',1e-6,@(x) validateattributes(x,{'numeric'},{'scalar'},'ModelOptL','TolCon'));
    p.addParameter('MaxIter',1000,@(x) validateattributes(x,{'numeric'},{'scalar','integer'},'ModelOptL','MaxIter'));
    p.addParameter('Algorithm','interior-point',@(x) parsevalidatestring(x,{'trust-region-reflective','active-set','interior-point','sqp'},'ModelOptL','Algorithm'));
    p.addParameter('Display','off',@(x) parsevalidatestring(x,{'off','iter','iter-detailed','notify','notify-detailed','final','final-detailed'},'ModelOptL','Display'));
    p.addParameter('Hessian','fcn',@(x) parsevalidatestring(x,{'fcn','mult','none'},'ModelOptL','Hessian'));
    p.addParameter('fp',0.5,@(x) validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1},'ModelOptL','fp'));
end
p.parse(varargin{:});
r=p.Results;
extra=[fields(p.Unmatched)'; struct2cell(p.Unmatched)'];

fp=r.fp;

n=length(Wp);
w=BinaryWeights(n);

x0 = Mats2Params(Wp,Wm);            % Starting guess 
lb=zeros(size(x0));
ub=ones(size(x0));
[linconstr_A,linconstr_b]=ParamsConstraints(n);

options = optimoptions('fmincon',extra{:},'Algorithm',r.Algorithm,'Display',r.Display,...
    'TolFun', r.TolFun,...  % termination based on function value (of the derivative)
    'TolX', r.TolX,...
    'TolCon',r.TolCon,...
    'MaxIter',r.MaxIter);

if r.UseDerivs
    options = optimoptions(options,'GradObj','on');
    if strcmpi(r.Hessian,'fcn')
        options=optimoptions('fmincon',options,'Hessian','user-supplied',...
            'HessFcn',@opthess);
    elseif strcmpi(r.Hessian,'mult')
        options=optimoptions('fmincon',options,'Hessian','user-supplied',...
            'SubproblemAlgorithm','cg','HessMult',@opthessmult);
    end
    [x,A,ef] = fmincon(@(y)OptFunGradL(y,sm,fp,w),x0,...
        linconstr_A,linconstr_b,...
        [],[],lb,ub,[],...
        options);%fun,xo,A,b,Aeq,beq,lb,ub,nonlcon,options
else
    [x,A,ef] = fmincon(@(y)OptFunL(y,sm,fp,w),x0,...
        linconstr_A,linconstr_b,...
        [],[],lb,ub,[],...
        options);%fun,xo,A,b,Aeq,beq,lb,ub,nonlcon,options
end
[Wp,Wm]=Params2Mats(x);

% [~,~,ix]=SortByEta(0.5*Wp+0.5*Wm,w);
[newWp,newWm]=SortByWtEtaS(Wp,Wm,w,fp,sm);
A=-A;

if r.DispExit
    disp(ExitFlagMsg(ef));
end


    function h=opthess(x,~)
        h = OptHessL(x,sm,fp,w);
    end

    function Hv=opthessmult(x,~,v)
        Hv = OptHessMultL(x,sm,fp,w,v);
    end


end

