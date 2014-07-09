function [ newobj ] = Build( func,varargin )
%newobj=BUILD(func,extraArgs,fp) build SynapseIdModel from function
%   [Wp,Wm,newobj.w]=func(extraArgs{:})
%   newobj.Initial=EqProb(fp*Wp+(1-fp)*Wm)
%   newobj.M={Wp+eye,Wm+eye}


persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='SynapseIdModel.Build';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addOptional('extraArgs',{});
    p.addOptional('fp',0.5);
    p.addParameter('ScaleW',1);
end
p.parse(varargin{:});
r=p.Results;

if ischar(func)
    func=str2func(func);
end

[Wp,Wm,w]=func(r.extraArgs{:});

newobj=SynapseIdModel;

newobj=newobj.setM({r.ScaleW*Wp+eye(length(Wp)),r.ScaleW*Wm+eye(length(Wp))});
newobj=newobj.setW(w);
newobj=newobj.setFp(r.fp);
newobj=newobj.CalcEqProb;


end

