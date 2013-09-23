function [ obj3 ] = mrdivide( obj1,obj2 )
%MRDIVIDE division of probs by scalar
%   only adds obj.M and obj.initial. other props taken from obj1

[obj3,thescalar]=extractArgOfType({obj1,obj2},'SynapseIdModel');
thescalar=thescalar{1};
obj3.initial=obj3.initial/thescalar;
for i=1:length(obj3.M)
    obj3.M{i}=obj3.M{i}/thescalar;
end


end

