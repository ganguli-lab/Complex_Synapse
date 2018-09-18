function [ Zinv,piv ] = GetZinv( obj,varargin )
%[Zinv,piv]=SynapseMemoryModel.GETZINV inverse of fundamental matrix
%   Zinv = ev*piv - Wf
%   ev = column vector of ones
%   piv = arbitrary row vec, s.t. piv*ev~=0 (default ev')

ev=ones(size(obj.w));

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='SynapseMemoryModel.GetZinv';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addOptional('piv',[],@(x)validateattributes(x,{'numeric'},{'row'},'SynapseMemoryModel.GetZinv','piv',2));
end
p.parse(varargin{:});
piv=p.Results.piv;
if isempty(piv)
    piv=ev';
end

Zinv = ev*piv - obj.GetWf;


end

