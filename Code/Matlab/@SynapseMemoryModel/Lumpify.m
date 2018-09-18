function [ newobj ] = Lumpify( obj,partitions )
%newobj=obj.LUMPIFY(partitions) lump states together
%   pertitions = cell array of vectors containing indices of states in each
%                partition

[U,V]=LumpProj(partitions);

w=U*obj.w;
Wp=U*obj.Wp*V;
Wm=U*obj.Wm*V;

newobj=obj.setWp(Wp);
newobj=newobj.setWm(Wm);
newobj=newobj.setW(w);

end

