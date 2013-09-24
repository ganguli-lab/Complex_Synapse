function [ tf ] = isvalid( obj )
%tf=ISVALID(obj) are all matrices of appropriate size and normalised?

tf=isprob(obj.Initial) &&...
    ~isempty(obj.M) &&...
    isstochasticD(obj.M{1}) &&...
    length(obj.Initial)==length(obj.M{1});
for i=2:length(obj.M)
    tf=tf && isstochasticD(obj.M{i}) &&...
    samesize(obj.M{i},obj.M{1});
end
for i=1:numel(obj.outProj)
    tf=tf && samesize(obj.outProj{i},obj.M{1});
end


end

