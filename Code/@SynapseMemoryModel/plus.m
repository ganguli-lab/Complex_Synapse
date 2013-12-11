function [ obj3 ] = plus( obj1,obj2 )
%PLUS addition of probs
%   only adds obj.M and obj.Initial. other props taken from obj1

error(CheckSize(obj2,@(x) obj1.SameSizes(x),'samesize(obj1)'));

obj3=obj1;
obj3.Wp=obj3.Wp+obj2.Wp;
obj3.Wm=obj3.Wm+obj2.Wm;

end

