function [ newobj,ix ] = Sort( obj )
%[newobj,ix]=SORT(obj) Put states in order of decreasing eta^+
%   obj,newobj = SynapseMemoryModel
%   ix = order of states, newobj.w=obj.w(ix)


[~,ix]=sort(obj.w,'ascend');
newobj=obj.Reorder(ix);

deta=obj.DeltaEta;

wchange=[0 find(diff(newobj.w)) length(newobj.w)];

ix=[];
for i=1:length(wchange)-1
    [~,dix]=sort(deta(wchange(i)+1:wchange(i+1)),'descend');
    ix=[ix dix'+wchange(i)];
end

newobj=obj.Reorder(ix);


end

