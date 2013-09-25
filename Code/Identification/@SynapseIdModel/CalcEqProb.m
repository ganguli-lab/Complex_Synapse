function newobj=CalcEqProb(obj,fp)
%newobj=CALCEQPROB(obj,fp) set Initial to eq dist

newobj=obj;
W=zeros(length(newobj.w));
fp(end+1)=1-sum(fp);

for i=1:length(obj.M)
    W=W+fp(i)*newobj.M{i};
end
newobj=newobj.setInitial(EqProb(W-eye(length(W))));

end

