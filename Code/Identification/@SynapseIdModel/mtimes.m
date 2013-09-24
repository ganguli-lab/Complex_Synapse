function [ obj3 ] = mtimes( obj1,obj2 )
%MTIMES multiplication of probs by scalar
%   only adds obj.M and obj.Initial. other props taken from obj1

[obj3,thescalar]=extractArgOfType({obj1,obj2},'SynapseIdModel');
thescalar=thescalar{1};
obj3.Initial=obj3.Initial*thescalar;
for i=1:length(obj3.M)
    obj3.M{i}=obj3.M{i}*thescalar;
end


end

