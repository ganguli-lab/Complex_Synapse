function [ newobj ] = Build( func,varargin )
%newobj=BUILD(func,fp,extraArgs) build SynapseMemoryModel from function
%   [newobj.Wp,newobj.Wm,newobj.w]=func(extraArgs{:})


persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='SynapseMemoryModel.Build';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addRequired('func',@(x)validateattributes(x,{'function_handle','char'},{},'SynapseMemoryModel.Build','func',1));
    p.addOptional('fp',0.5,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1},'SynapseMemoryModel.Build','fp',2));
    p.addOptional('extraArgs',{},@(x)validateattributes(x,{'cell'},{},'SynapseMemoryModel.Build','extraArgs',3));
    p.addParameter('ScaleW',1,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'},'SynapseMemoryModel.Build','ScaleW'));
end
p.parse(func,varargin{:});
r=p.Results;



if ischar(r.func)
    func=str2func(r.func);
end

[Wp,Wm,w]=func(r.extraArgs{:});

newobj=SynapseMemoryModel;

newobj=newobj.setWp(r.ScaleW*Wp);
newobj=newobj.setWm(r.ScaleW*Wm);
newobj=newobj.setFp(r.fp);
newobj=newobj.setW(w);

end

