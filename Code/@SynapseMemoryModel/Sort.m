function [ newobj,ix ] = Sort( obj )
%newobj=SORT(obj) Put states in order of decreasing eta^+


newobj=obj;

[newobj.w,ix]=sort(obj.w,'ascend');
newobj.Wp=obj.Wp(ix,ix);
newobj.Wm=obj.Wm(ix,ix);

deta=obj.DeltaEta;

wchange=[0 find(diff(newobj.w)) length(newobj.w)];

ix=[];
for i=1:length(wchange)-1
    [~,dix]=sort(deta(wchange(i)+1:wchange(i+1)),'descend');
    ix=[ix dix'+wchange(i)];
end

newobj.Wp=obj.Wp(ix,ix);
newobj.Wm=obj.Wm(ix,ix);


end

