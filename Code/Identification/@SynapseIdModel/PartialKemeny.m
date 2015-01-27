function [ deta ] = PartialKemeny( obj )
%deta=obj.PARTIALKEMENY partial kemeny "constants"
%   deta_i = eta^+_i - eta^+_n

W=zeros(obj.NumStates);
fp=obj.fp;
fp(end+1)=1-sum(fp);

for i=1:length(obj.M)
    W=W+fp(i)*obj.M{i};
end

deta=DeltaEta(W-eye(obj.NumStates),obj.w);

end

