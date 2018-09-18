function [ tf ] = MatchW( obj,otherobj )
%tf=MATCHW(obj,otherobj) do obj and other obj have the same w?

tf = obj.NumStates==otherobj.NumStates && all(obj.w == otherobj.w);

end

