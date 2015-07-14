function [ orphans ] = FindOrphans( obj )
%orphans=obj.FINDORPHANS find states that are never visited
%   orphans = row vector of never visited state indices
%   obj = SynapseMemoryModel


Wf=obj.GetWf;

orphans=false(1,obj.NumStates);

for i=1:obj.NumStates
    inprob=sum(Wf(~orphans,:),1)-diag(Wf)';
    neworphans=(inprob<obj.OrphanThresh);
    if any(neworphans & ~orphans)
        orphans=orphans | neworphans;
    else
        break;
    end
end

orphans=find(orphans);

end

