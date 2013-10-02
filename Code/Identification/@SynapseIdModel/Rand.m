function [ newobj ] = Rand( w,varargin )
%newobj=Rand(obj,func,fp,varargin) build SynapseIdModel from function
%   [Wp,Wm,newobj.w]=func(varargin{:})
%   newobj.Initial=EqProb(fp*Wp+(1-fp)*Wm)
%   newobj.M={Wp+eye,Wm+eye}


ScaleW=1;
sparsity=1;
NumPlastTypes=2;
varargin=assignApplicable(varargin);

M=cell(1,NumPlastTypes);
for i=1:NumPlastTypes
    M{i}=ScaleW*RandTrans(length(w),sparsity)+eye(length(w));
end
p=rand(1,length(w));


newobj=SynapseIdModel;

newobj=newobj.setM(M);
newobj=newobj.setW(w);
newobj=newobj.setInitial(p/sum(p));


end

