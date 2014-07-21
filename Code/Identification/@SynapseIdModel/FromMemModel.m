function [ idmodelobj ] = FromMemModel( memmodelobj,varargin )
%idmodelobj=SynapseMemoryModel.FROMIDMODEL(memmodelobj,oldidmodelobj) build
%SynapseIdModel from SynapseMemoryModel
%   memmodelobj.Wp = idmodelobj.M{1}-eye
%   memmodelobj.Wm = idmodelobj.M{2}-eye
%   memmodelobj.w  = idmodelobj.w
%   memmodelobj.fp = idmodelobj.fp
%   optional arg: oldidmodelobj, if you want to copy other props

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='SynapseIdModel.FromMemModel';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addOptional('oldidmodelobj',SynapseIdModel,@(x)validateattributes(x,{'SynapseIdModel'},{},'SynapseIdModel.FromMemModel','oldidmodelobj',2));
end
p.parse(varargin{:});


idmodelobj=p.Results.oldidmodelobj;

idmodelobj=idmodelobj.setM({memmodelobj.Wp+eye(memmodelobj.NumStates),memmodelobj.Wm+eye(memmodelobj.NumStates)});
idmodelobj=idmodelobj.setW(memmodelobj.w);
idmodelobj=idmodelobj.setFp(memmodelobj.fp);
idmodelobj=idmodelobj.CalcEqProb;
idmodelobj=idmodelobj.SetWValInds;

end

