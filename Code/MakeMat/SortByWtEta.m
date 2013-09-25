function [ newM,neww,ix ] = SortByWtEta( M,w,fp )
%[newW,neww,ix]=SORTBYETA(M,w,fp) Put states in order of decreasing eta^+
%   M = transition rates, {Mpot,Mdep}
%   w = Weights of states (+/-1)
%   ix=sort order
%   newW=W(ix,ix)
%   neww=w(ix)

error(CheckType(M,'cell'));
error(CheckSize(M,@(x) numel(x)==2,'2 elements'));
error(CheckSize(M{1},@isstochasticD));
error(CheckSize(M{2},@isstochasticD));
error(CheckSize(M{2},@(x) samesize(x,M{1}),'samesize(Mpot)'));
error(CheckSize(w,@(x) length(x)==length(M{1}),'samesize(M)'));

[neww,ix]=sort(w,'ascend');
newM={M{1}(ix,ix),M{2}(ix,ix)};

deta=DeltaEta(fp*newM{1}+(1-fp)*newM{2}-eye(length(w)),w);

wchange=[0 find(diff(w)) length(w)];

ix=[];
for i=1:length(wchange)-1
    [~,dix]=sort(deta(wchange(i)+1:wchange(i+1)),'descend');
    ix=[ix dix'+wchange(i)];
end

newM={newM{1}(ix,ix),newM{2}(ix,ix)};
neww=neww(ix);


end