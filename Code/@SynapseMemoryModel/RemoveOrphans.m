function [ newobj ] = RemoveOrphans( obj,orphans )
%newobj=obj.REMOVEORPHANS(orphans) remove states that are never visited
%   orphans = row vector of never visited state indices
%   obj,newobj = SynapseMemoryModel


newobj=obj;

newobj.Wp(orphans,:)=[];
newobj.Wp(:,orphans)=[];

newobj.Wm(orphans,:)=[];
newobj.Wm(:,orphans)=[];

newobj.w(orphans)=[];

newobj=newobj.Normalise;

end

