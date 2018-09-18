function tf=SameSizes(obj,otherobj)
%tf=SAMESIZES(obj,otherobj) do the properties of obj and otherobj have the
%same sizes?
%   NumT

tf = obj.NumT == otherobj.NumT;


end

