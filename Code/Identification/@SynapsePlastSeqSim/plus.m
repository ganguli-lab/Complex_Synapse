function [ obj3 ] = plus( obj1,obj2 )
%PLUS concatenation of properties
%   concatenates obj.potdep, obj.stateseq and obj.readouts. other props taken from obj1

obj3=plus@SynapsePlastSeq(obj1,obj2);
obj3.stateseq=[obj3.stateseq obj2.stateseq];


end

