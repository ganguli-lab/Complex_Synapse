function [ obj3 ] = plus( obj1,obj2 )
%PLUS addition of probs
%   only adds obj.M and obj.Initial. other props taken from obj1

% error(CheckSize(obj2,@(x) obj1.SameSizes(x),'samesize(obj1)'));

obj3=obj1;
obj3.Initial=obj3.Initial+obj2.Initial;
for i=1:length(obj3.M)
    obj3.M{i}=obj3.M{i}+obj2.M{i};
end


end

