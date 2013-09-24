function [ newobj ] = Zero( obj )
%ZERO set all probs to zero
%   M{:} and initial

newobj=obj;
newobj.initial=zeros(size(newobj.initial));
for i=1:length(newobj.M)
    newobj.M{i}=zeros(size(newobj.M{i}));
end

end

