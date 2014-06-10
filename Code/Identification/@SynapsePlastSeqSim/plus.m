function [ obj3 ] = plus( obj1,obj2 )
%PLUS concatenation of properties
%   concatenates obj.potdep, obj.stateseq and obj.readouts. other props taken from obj1

obj3=obj1;
obj3.potdep=[obj3.potdep obj2.potdep];
obj3.stateseq=[obj3.stateseq obj2.stateseq];
obj3.readouts=[obj3.readouts obj2.readouts];


end

