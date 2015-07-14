function [ newobj ] = Reorder( obj,ix )
%newobj=obj.REORDER(ix) change order of states
%   obj,newobj = SynapseMemoryModel
%   ix = order of states, newobj.w=obj.w(ix)

newobj=obj;

newobj.Wp=obj.Wp(ix,ix);
newobj.Wm=obj.Wm(ix,ix);

newobj.w=obj.w(ix);

end

