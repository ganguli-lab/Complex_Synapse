function [ newobj ] = Build( func,fp,varargin )
%newobj=BUILD(func,fp,varargin) build SynapseIdModel from function
%   [Wp,Wm,newobj.w]=func(varargin{:})
%   newobj.Initial=EqProb(fp*Wp+(1-fp)*Wm)
%   newobj.M={Wp+eye,Wm+eye}


ScaleW=1;
varargin=assignApplicable(varargin);

if ischar(func)
    func=str2func(func);
end

[Wp,Wm,w]=func(varargin{:});

newobj=SynapseIdModel;

newobj=newobj.setM({ScaleW*Wp+eye(length(Wp)),ScaleW*Wm+eye(length(Wp))});
newobj=newobj.setW(w);
newobj=newobj.setFp(fp);
newobj=newobj.CalcEqProb;


end

