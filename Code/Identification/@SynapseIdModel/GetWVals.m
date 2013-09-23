function [ wvals ] = GetWVals( obj )
%wvals=GETWVALS(obj) list of different values in obj.w

wvals=sort(obj.w);
wvals([diff(wvals)==0;false])=[];


end

