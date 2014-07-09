function [ newobj ] = Rand( w,fp,varargin )
%newobj=Rand(obj,func,fp,varargin) build SynapseIdModel from function
%   [Wp,Wm,newobj.w]=func(varargin{:})
%   newobj.Initial=EqProb(fp*Wp+(1-fp)*Wm)
%   newobj.M={Wp+eye,Wm+eye}


ScaleW=1;
sparsity=1;
varargin=assignApplicable(varargin);


newobj=SynapseIdModel;

newobj=newobj.setWp(ScaleW*RandTrans(length(w),sparsity));
newobj=newobj.setWm(ScaleW*RandTrans(length(w),sparsity));
newobj=newobj.setFp(fp);
newobj=newobj.setW(w);


end

