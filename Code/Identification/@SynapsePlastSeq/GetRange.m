function newobj=GetRange(obj,range)
%newobj=GETRANGE(obj,range) extract subset of SynapsePlastSeq for
%range(1)<=t<=range(2)

newobj=obj;

newobj.potdep=newobj.potdep(range(1):range(2));
newobj.stateseq=newobj.stateseq(range(1):range(2));
newobj.readouts=newobj.readouts(range(1):range(2));


end

