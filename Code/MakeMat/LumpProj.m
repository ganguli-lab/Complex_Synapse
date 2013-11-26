function [ U,V ] = LumpProj( partitions )
%[U,V]=LUMPPROJ(num_states,partitions) matrices for testing lumpability and
%constructing lumped process
%   U: Ux = average x over subsets
%   V: pV = sum of p over subsets
%   partitions = cell array of vectors containing indices of states in each
%                partition

error(CheckType(partitions,'cell'));

V=zeros(max([partitions{:}]),length(partitions));

for a=1:length(partitions)
    for j=1:length(partitions{a})
        V(partitions{a}(j),a)=1;
    end
end

U=diag(1./sum(V,1))*V';


end

