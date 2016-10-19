function [ tf ] = isvalid( obj )
%tf=ISVALID(obj) are all matrices of appropriate size and normalised?

tf= isvalid@VORtrainSeq(obj) &&...
    isrow(obj.fps_other) && length(obj.fps_other)==length(obj.tTrain)+1;


end

