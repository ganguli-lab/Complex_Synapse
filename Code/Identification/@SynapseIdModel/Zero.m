function [ newobj ] = Zero( obj )
%ZERO set all probs to zero
%   M{:} and Initial

newobj=obj;
newobj.Initial=zeros(size(newobj.Initial));
for i=1:length(newobj.M)
    newobj.M{i}=zeros(size(newobj.M{i}));
end

end

