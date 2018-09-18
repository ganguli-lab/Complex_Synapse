function tf=SameSizes(obj,otherobj)
%tf=SAMESIZES(obj,otherobj) do the properties of obj and otherobj have the
%same sizes?
%   NUmStates, NumPlast and NumWvals

tf = obj.NumStates == otherobj.NumStates &&...
    obj.NumPlast == otherobj.NumPlast &&...
    obj.NumWvals == otherobj.NumWvals;


end

