function [ tf ] = isvalid( obj )
%tf=ISVALID(obj) are all matrices of appropriate size and normalised?

tf = obj.VORrel.isvalid && obj.VORcomp.isvalid &&...
    all(obj.VORrel.tTrain == obj.VORcomp.tTrain);


end

