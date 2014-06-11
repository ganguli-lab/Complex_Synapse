function newobj=GetRange(obj,range)
%newobj=GETRANGE(obj,range) extract subset of SynapsePlastSeq for
%range(1)<=t<=range(2)

newobj=GetRange@SynapsePlastSeq(obj,range);

newobj.stateseq=newobj.stateseq(range(1):range(2));


end

