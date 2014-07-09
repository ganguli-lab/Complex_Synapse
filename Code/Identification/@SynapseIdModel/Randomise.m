function [ newobj ] = Randomise( obj,varargin )
%newobj=obj.RANDOMISE(w,...) set Initial and M to random
%   w = new synaptic weights (optional)
%   parameter value pairs (default):
%   'NumPlastTypes' number of typse of plasticity (2),
%   'ScaleW' multiplies elelments (1),
%   'sparsity' probability that each elemnt is non-zero (1).

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='SynapseIdModel.Randomise';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addOptional('w',obj.w);
    p.addParameter('NumPlastTypes',obj.NumPlast);
    p.addParameter('ScaleW',1);
    p.addParameter('sparsity',1);
%     p.addParameter('Normalise',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
end
p.parse(varargin{:});
if any(strcmp('w',p.UsingDefaults))
    w=obj.w;
else
    w=p.Results.w;
end

newobj=obj.setW(w);

init=rand(1,length(newobj.w));
newobj=newobj.setInitial(init/sum(init));

M=cell(1,p.Results.NumPlastTypes);
for i=1:length(M)
    M{i}=p.Results.ScaleW*RandTrans(newobj.NumStates,p.Results.sparsity)+eye(newobj.NumStates);
end
newobj=newobj.setM(M);

end

