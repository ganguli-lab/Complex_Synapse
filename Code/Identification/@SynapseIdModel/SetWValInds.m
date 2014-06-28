function [ newobj ] = SetWValInds( obj )
%newobj=obj.SETWVALINDS replace values of obj.w with ind of obj.outProj
%   Detailed explanation goes here

newobj=obj.setW(obj.GetWValInds);

end

