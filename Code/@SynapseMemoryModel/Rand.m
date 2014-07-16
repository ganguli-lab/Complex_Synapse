function [ newobj ] = Rand( w,varargin )
%newobj=Rand(obj,func,fp,varargin) build SynapseMemoryModel from function
%   [Wp,Wm,newobj.w]=func(varargin{:})
%   newobj.Initial=EqProb(fp*Wp+(1-fp)*Wm)
%   newobj.M={Wp+eye,Wm+eye}

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='SynapseMemoryModel.Rand';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addRequired('w',@(x)validateattributes(x,{'numeric'},{'column'},'SynapseMemoryModel.Rand','w',1));
    p.addOptional('fp',0.5,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1},'SynapseMemoryModel.Rand','fp',2));
    p.addParameter('ScaleW',{'pot','dep'},@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'},'SynapseMemoryModel.Rand','ScaleW'));
    p.addParameter('sparsity',1,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1},'SynapseMemoryModel.Rand','sparsity'));
end
p.parse(w,varargin{:});
r=p.Results;



newobj=SynapseMemoryModel(p.Unmatched);

newobj=newobj.setWp(r.ScaleW*RandTrans(length(w),r.sparsity));
newobj=newobj.setWm(r.ScaleW*RandTrans(length(w),r.sparsity));
newobj=newobj.setFp(r.fp);
newobj=newobj.setW(r.w);


end

