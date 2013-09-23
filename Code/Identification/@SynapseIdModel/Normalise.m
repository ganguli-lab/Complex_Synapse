function [ newobj ] = Normalise( obj )
%newobj=NORMALISE(obj) normalise probs in obj (M & initial)

newobj=obj;

newobj.initial=newobj.initial/sum(newobj.initial);
for i=1:length(newobj.M)
    newobj.M{i}=StochastifyD(newobj.M{i});
end

end

