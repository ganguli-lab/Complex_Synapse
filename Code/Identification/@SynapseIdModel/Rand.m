function [ newobj ] = Rand( w,varargin )
%newobj=Rand(obj,func,fp,varargin) build SynapseIdModel from function
%   [Wp,Wm,newobj.w]=func(varargin{:})
%   newobj.Initial=EqProb(fp*Wp+(1-fp)*Wm)
%   newobj.M={Wp+eye,Wm+eye}


ScaleW=1;
sparsity=1;
varargin=assignApplicable(varargin);

Wp=RandTrans(length(w),sparsity);
Wm=RandTrans(length(w),sparsity);
p=rand(1,length(w));


newobj=SynapseIdModel;

newobj=newobj.setM({ScaleW*Wp+eye(length(Wp)),ScaleW*Wm+eye(length(Wp))});
newobj=newobj.setW(w);
newobj=newobj.setInitial(p/sum(p));


end

