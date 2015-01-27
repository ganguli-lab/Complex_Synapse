function [ partitions ] = FindLumps( obj )
%partitions=obj.FINDLUMPS find partitions of states that could be lumpable
%   pertitions = cell array of vectors containing indices of states in each
%                partition
%   checks if eta^\pm are degenerate to within obj.DegThresh.

[newobj,ix]=obj.Sort;


deta=newobj.PartialKemeny;

detachange=[0; find( abs(diff(deta))>newobj.DegThresh | diff(newobj.w)); newobj.NumStates];
detachange=diff(detachange);


partitions=mat2cell(ix,1,detachange);

end

