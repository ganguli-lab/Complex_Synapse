function [ newWp,newWm,A,ef ] = ModelOptDoubleIL( Wp,Wm,sm,sc,Ac,varargin)
%[newWp,newWm,A,ef]=MODELOPTDOUBLEIL(Wp,Wm,sm) run gradient descent on model to
%maximise A(sm) subject to A(sc)>=Ac
%   sm = inverse time value
%   WP = potentiation transition rates
%   WM = depression transition rates
%   A  = Laplace Transf value
%   ef = exit flag

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='ModelOptDoubleIL';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addParameter('UseDerivs',true,@(x) validateattributes(x,{'logical'},{'scalar'},'ModelOptDoubleIL','UseDerivs'));
    p.addParameter('DispExit',false,@(x) validateattributes(x,{'logical'},{'scalar'},'ModelOptDoubleIL','DispExit'));
    p.addParameter('TolFun',1e-6,@(x) validateattributes(x,{'numeric'},{'scalar'},'ModelOptDoubleIL','TolFun'));
    p.addParameter('TolX',1e-10,@(x) validateattributes(x,{'numeric'},{'scalar'},'ModelOptDoubleIL','TolFun'));
    p.addParameter('TolCon',1e-6,@(x) validateattributes(x,{'numeric'},{'scalar'},'ModelOptDoubleIL','TolFun'));
    p.addParameter('MaxIter',1000,@(x) validateattributes(x,{'numeric'},{'scalar','integer'},'ModelOptDoubleIL','TolFun'));
    p.addParameter('Algorithm','interior-point',@(x) validatestring(x,{'trust-region-reflective','active-set','interior-point','sqp'},'ModelOptDoubleIL','TolFun'));
    p.addParameter('Display','off',@(x) validatestring(x,{'off','iter','iter-detailed','notify','notify-detailed','final','final-detailed'},'ModelOptDoubleIL','TolFun'));
    p.addParameter('fp',0.5,@(x) validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1},'ModelOptDoubleIL','fp'));
end
p.parse(varargin{:});
r=p.Results;

fp=r.fp;

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
    options = optimset(options,'GradObj','on',...
        'Hessian','user-supplied','HessFcn',@lagrangianhessmax);
    [x1] = fmincon(@initoptfngr,x0,...
        linconstr_A,linconstr_b,...
        [],[],lb,ub,...
        [],... 
        options);%fun,xo,A,b,Aeq,beq,lb,ub,nonlcon,options
else
    [x1] = fmincon(@initoptfn,x0,...
        linconstr_A,linconstr_b,...
        [],[],lb,ub,...
        [],... 
        options);%fun,xo,A,b,Aeq,beq,lb,ub,nonlcon,options
end
% x1=x0;

if r.UseDerivs
    options = optimset(options,'GradObj','on','GradConstr','on',...
        'Hessian','user-supplied','HessFcn',@lagrangianhess);
    [x,A,ef] = fmincon(@(y)OptFunGradL(y,sm,fp,w),x1,...
        linconstr_A,linconstr_b,...
        [],[],lb,ub,...
        @nlconstrgr,... 
        options);%fun,xo,A,b,Aeq,beq,lb,ub,nonlcon,options
else
    [x,A,ef] = fmincon(@(y)OptFunL(y,sm,fp,w),x1,...
        linconstr_A,linconstr_b,...
        [],[],lb,ub,...
        @nlconstr,... 
        options);%fun,xo,A,b,Aeq,beq,lb,ub,nonlcon,options
end

% if any(linconstr_A*x>linconstr_b)
%     A=0;
% end
    

[Wp,Wm]=Params2Mats(x);

% [~,~,ix]=SortByEta(0.5*Wp+0.5*Wm,w);
[newWp,newWm]=SortByWtEtaS(Wp,Wm,w,fp,sm);
A=-A;

% [~,~,~,~,fail]=DoubleLaplace(sm,newWp,newWm,fp,w);

if ef==-2 || nlconstr(x)>0
    A=0;
end

if r.DispExit
    disp(ExitFlagMsg(ef));
end



    function [c,ceq]=nlconstr(xc)
        ceq=[];
        c=OptFunL(xc,sc,fp,w)+Ac;
    end

    function [c,ceq,gradc,gradceq]=nlconstrgr(xc)
        ceq=[];
        [c,gradc]=OptFunGradL(xc,sc,fp,w);
        c=c+Ac;
        gradceq=[];
    end

    function f=initoptfn(xi)
        f=OptFunL(xi,sc,fp,w);
%         f=OptFunL(xi,sc,fp,w)+Ac;
%         f=0.5*f^2;
    end

    function [f,gr]=initoptfngr(xi)
        [f,gr]=OptFunGradL(xi,sc,fp,w);
%         f=f+Ac;
%         gr=f*gr;
%         f=0.5*f^2;
    end

    function h=lagrangianhess(x,lambda)
        h = OptHessL(x,sm,fp,w) + lambda.ineqnonlin * OptHessL(x,sc,fp,w);
    end

    function h=lagrangianhessmax(x,~)
        h = OptHessL(x,sm,fp,w);
    end

end

