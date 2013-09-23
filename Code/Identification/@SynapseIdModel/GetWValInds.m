function [ wvalinds ] = GetWValInds( obj )
%wvals=GETWVALINDS(obj) obj.w with values replaced by ind of obj.outProj

wvals=obj.GetWVals;
wvalinds=obj.w;

for i=1:length(wvals);
    wvalinds(obj.w==wvals(i))=i;
end


end

