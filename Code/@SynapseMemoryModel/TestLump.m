function [ tf,Mtest,wtest ] = TestLump( obj,partitions,varargin )
%tf=obj.TESTLUMP(partitions) test if SynapseIdModel obj is lumpabel wrt
%partions
%   pertitions = cell array of vectors containing indices of states in each
%                partition
%   checks if VUw=w and VUMV=MV to within obj.LumpThresh
%   V,U from LumpProj. Uses LumpTest.

% persistent p
% if isempty(p)
%     p=inputParser;
%     p.FunctionName='SynapseMemoryModel.Spectrum';
%     p.StructExpand=true;
%     p.KeepUnmatched=false;
%     p.addParameter('LumpThresh',1e-3,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'},'SynapseMemoryModel.TestLump','LumpThresh'));
% end
% p.parse(varargin{:});
% r=p.Results;


wtest=LumpTest(obj.w,partitions);
wtest=max(abs(wtest));
Mtest=cellfun(@(x)LumpTest(x,partitions),{obj.Wp,obj.Wm},'UniformOutput',false);
Mtest=cellfun(@(x)max(max(abs(x))),Mtest);


tf=max([wtest Mtest])<obj.LumpThresh;


end

