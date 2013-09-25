function [ newobj,ix ] = Sort( obj,fp )
%newobj=SORT(obj) Put states in order of decreasing eta^+

error(CheckSize(fp,@(x) length(x)==length(obj.M)-1,'length(M)-1'));

newobj=obj;

W=zeros(length(newobj.w));
fp(end+1)=1-sum(fp);

[newobj.w,ix]=sort(obj.w,'ascend');
newobj.Initial=newobj.Initial(ix);
for i=1:length(obj.M)
    newobj.M{i}=obj.M{i}(ix,ix);
    W=W+fp(i)*newobj.M{i};
end

deta=DeltaEta(W-eye(length(newobj.w)),newobj.w);

wchange=[0 find(diff(newobj.w)) length(newobj.w)];

ix=[];
for i=1:length(wchange)-1
    [~,dix]=sort(deta(wchange(i)+1:wchange(i+1)),'descend');
    ix=[ix dix'+wchange(i)];
end

for i=1:length(obj.M)
    newobj.M{i}=newobj.M{i}(ix,ix);
end
newobj.Initial=newobj.Initial(ix);
newobj=newobj.CalcOutProj;


end

