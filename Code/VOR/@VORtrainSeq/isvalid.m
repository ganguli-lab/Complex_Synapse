function [ tf ] = isvalid( obj )
%tf=ISVALID(obj) are all matrices of appropriate size and normalised?

tf = isrow(obj.tTrain) &&...
    isrow(obj.fps) && length(obj.fps)==length(obj.tTrain)+1 &&...
    isrow(obj.rs) && length(obj.rs)==length(obj.tTrain)+1;


end

