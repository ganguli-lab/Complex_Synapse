function [ obj3 ] = minus( obj1,obj2 )
%PLUS subtraction of probs
%   only subtracts obj.M and obj.Wp. other props taken from obj1

error(CheckSize(obj2,@(x) obj1.SameSizes(x),'samesize(obj1)'));

obj3=obj1;
obj3.Wp=obj3.Wp-obj2.Wp;
obj3.Wm=obj3.Wm-obj2.Wm;


end

