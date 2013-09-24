function [ newobj ] = Build( func,fp,varargin )
%newobj=BUILD(obj,func,fp,varargin) build SynapseIdModel from function
%   [Wp,Wm,newobj.w]=func(varargin{:})
%   newobj.Initial=EqProb(fp*Wp+(1-fp)*Wm)
%   newobj.M={Wp+eye,Wm+eye}

if ischar(func)
    func=str2func(func);
end

[Wp,Wm,w]=func(varargin{:});

newobj=SynapseIdModel;

newobj=newobj.setInitial(EqProb(fp*Wp+(1-fp)*Wm));
newobj=newobj.setM({Wp+eye(length(Wp)),Wm+eye(length(Wp))});
newobj=newobj.setW(w);


end

