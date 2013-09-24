function [ newobj ] = Normalise( obj )
%newobj=NORMALISE(obj) normalise probs in obj (M & Initial)

newobj=obj;

newobj.Initial=newobj.Initial/sum(newobj.Initial);
for i=1:length(newobj.M)
    newobj.M{i}=StochastifyD(newobj.M{i});
end

end

