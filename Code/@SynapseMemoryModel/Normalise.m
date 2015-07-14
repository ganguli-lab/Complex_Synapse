function [ newobj ] = Normalise( obj )
%newobj=obj.NORMALISE normalise transition probabilites
%   ensure that row sums are zero
%   obj,newobj = SynapseMemoryModel

newobj=obj;
newobj.Wp=StochastifyC(newobj.Wp);
newobj.Wm=StochastifyC(newobj.Wm);

end

