function [ partitions ] = FindLumps( obj,varargin )
%partitions=obj.FINDLUMPS find partitions of states that could be lumpable
%   pertitions = cell array of vectors containing indices of states in each
%                partition
%   checks if eta^\pm are degenerate to within obj.DegThresh.

% persistent p
% if isempty(p)
%     p=inputParser;
%     p.FunctionName='SynapseMemoryModel.Spectrum';
%     p.StructExpand=true;
%     p.KeepUnmatched=false;
%     p.addParameter('DegThresh',1e-3,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'},'SynapseMemoryModel.FindLumps','DegThresh'));
% end
% p.parse(varargin{:});
% r=p.Results;


[newobj,ix]=obj.Sort;


deta=newobj.DeltaEta;

detachange=[0; find( abs(diff(deta))>obj.DegThresh | diff(newobj.w)); newobj.NumStates];
detachange=diff(detachange);


partitions=mat2cell(ix,1,detachange);

end

