function [ partitions ] = Str2partitions( part_str )
%partitions=STR2PARTITIONS(part_str) construct partitions from string
%   partitions = cell array of vectors containing indices of states in each
%                partition
%   part_str = lisiting states with partitions separated by semi-colons

error(CheckType(part_str,'char'));

if part_str(end)~=';'
    part_str=[part_str ';'];
end

part_bnd = find(part_str==';');

partitions = cell(1,length(part_bnd));

part_bnd = [0 part_bnd];

for i=1:length(partitions)
    partitions{i} = eval(['[' part_str((part_bnd(i)+1):part_bnd(i+1)) ']']);
end


end

