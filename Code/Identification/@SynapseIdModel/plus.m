function [ obj3 ] = plus( obj1,obj2 )
%PLUS addition of probs
%   only adds obj.M and obj.initial. other props taken from obj1

error(CheckSize(obj2,@(x) samesize(x.initial,obj1.initial),'samesize(obj1)'));
error(CheckSize(obj2,@(x) samesize(x.M,obj1.M),'samesize(obj1)'));
obj3=obj1;
obj3.initial=obj3.initial+obj2.initial;
for i=1:length(obj3.M)
    obj3.M{i}=obj3.M{i}+obj2.M{i};
end


end

