function [ tf ] = TestLump( obj,partitions )
%tf=obj.TESTLUMP(partitions) test if SynapseIdModel obj is lumpabel wrt
%partions
%   pertitions = cell array of vectors containing indices of states in each
%                partition
%   checks if VUw=w and VUMV=MV to within obj.LumpThresh
%   V,U from LumpProj. Uses LumpTest.


wtest=LumpTest(obj.w,partitions);
wtest=max(abs(wtest));
Mtest=cellfun(@(x)LumpTest(x,partitions),obj.M,'UniformOutput',false);
Mtest=cellfun(@(x)max(max(abs(x))),Mtest);


tf=max([wtest Mtest])<obj.LumpThresh;


end

