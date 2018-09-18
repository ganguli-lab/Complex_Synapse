function [ newobj ] = Rand( w,varargin )
%newobj=Rand(w,...) build random SynapseIdModel
%   w = synaptic weights
%   parameter value pairs (default):
%   'NumPlastTypes' number of typse of plasticity (2),
%   'ScaleW' multiplies elelments (1),
%   'sparsity' probability that each elemnt is non-zero (1).
%other params passed to SynapseIdModel constructor


persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='SynapseIdModel.Rand';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addParameter('ScaleW',1);
    p.addParameter('sparsity',1);
    p.addParameter('NumPlastTypes',2);
%     p.addParameter('Normalise',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
end
p.parse(varargin{:});
r=p.Results;

M=cell(1,r.NumPlastTypes);
for i=1:r.NumPlastTypes
    M{i}=r.ScaleW*RandTrans(length(w),r.sparsity)+eye(length(w));
end
init=rand(1,length(w));


newobj=SynapseIdModel(p.Unmatched);

newobj=newobj.setM(M);
newobj=newobj.setW(w);
newobj=newobj.setInitial(init/sum(init));


end

