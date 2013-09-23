function [ tf ] = isvalid( obj )
%tf=ISVALID(obj) are all matrices of appropriate size and integers?

tf=isvector(obj.potdep) &&...
    isvector(obj.stateseq) &&...
    isvector(obj.readouts) &&...
    all(isint(obj.potdep)) &&...
    all(isint(obj.stateseq)) &&...
    all(isint(obj.readouts)) &&...
    length(obj.potdep)+1==length(obj.stateseq) &&...
    length(obj.stateseq)==length(obj.readouts);



end

