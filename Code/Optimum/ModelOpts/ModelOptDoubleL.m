function [ newWp,newWm,A,ef ] = ModelOptDoubleL( Wp,Wm,sm,sc,Ac,varargin)
%[newWp,newWm,A,ef]=MODELOPTDOUBLEL(Wp,Wm,sm) run gradient descent on model to
%maximise A(sm) subject to A(sc)=Ac
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

if r.UseDerivs
    options = optimset(options,'GradObj','on','GradConstr','on');
    [x,A,ef] = fmincon(@(y)OptFunGradDoubleL(y,sm,fp,w),x1,...
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

[~,~,~,~,fail]=DoubleLaplace(sm,newWp,newWm,fp,w);

if fail || ef==-2
    A=0;
end

if r.DispExit
    disp(ExitFlagMsg(ef));
end



    function [c,ceq]=nlconstr(xc)
        c=-1;
        ceq=OptFunL(xc,sc,fp,w)+Ac;
    end

    function [c,ceq,gradc,gradceq]=nlconstrgr(xc)
        c=-1;
        [ceq,gradceq]=OptFunGradL(xc,sc,fp,w);
        ceq=ceq+Ac;
        gradc=zeros(size(gradceq));
    end

    function f=initoptfn(xi)
        [~,f]=nlconstr(xi);
        f=0.5*f^2;
    end

    function [f,gr]=initoptfngr(xi)
        [~,f,~,gr]=nlconstrgr(xi);
        f=0.5*f^2;
        gr=f*gr;
    end

end

