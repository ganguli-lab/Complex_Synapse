function [ newobj ] = Rand( w,varargin )
%newobj=Rand(obj,func,fp,varargin) build SynapseIdModel from function
%   [Wp,Wm,newobj.w]=func(varargin{:})
%   newobj.Initial=EqProb(fp*Wp+(1-fp)*Wm)
%   newobj.M={Wp+eye,Wm+eye}


% ScaleW=1;
% sparsity=1;
% NumPlastTypes=2;
% varargin=assignApplicable(varargin);
persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='BWupdate';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addParameter('ScaleW',1);
    p.addParameter('sparsity',1);
    p.addParameter('NumPlastTypes',2);
%     p.addParameter('Normalise',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
end
p.parse(varargin{:});

M=cell(1,p.Results.NumPlastTypes);
for i=1:p.Results.NumPlastTypes
    M{i}=p.Results.ScaleW*RandTrans(length(w),p.Results.sparsity)+eye(length(w));
end
init=rand(1,length(w));


newobj=SynapseIdModel;

newobj=newobj.setM(M);
newobj=newobj.setW(w);
newobj=newobj.setInitial(init/sum(init));


end

