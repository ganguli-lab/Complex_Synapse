function [ tf ] = isvalid( obj )
%tf=ISVALID(obj) are all matrices of appropriate size and normalised?

tf= isstochasticC(obj.Wp) &&...
    isstochasticC(obj.Wm) &&...
    samesize(Wp,Wm) &&...
    length(obj.w)==length(obj.Wp) &&...
    inrange(obj.fp,0,1);


end

