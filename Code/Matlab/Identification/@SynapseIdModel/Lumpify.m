function [ newobj ] = Lumpify( obj,partitions )
%newobj=obj.LUMPIFY(partitions) lump states together
%   pertitions = cell array of vectors containing indices of states in each
%                partition

[U,V]=LumpProj(partitions);

w=U*obj.w;
Initial=obj.Initial*V;
M=cellfun(@(x) U*x*V, obj.M, 'UniformOutput',false);

newobj=obj.setInitial(Initial);
newobj=newobj.setM(M);
newobj=newobj.setW(w);

end

