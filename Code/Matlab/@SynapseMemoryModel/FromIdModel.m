function [ memmodelobj ] = FromIdModel( idmodelobj )
%memmodelobj=SynapseMemoryModel.FROMIDMODEL(idmodelobj) build
%SynapseMemoryModel from SynapseIdModel
%   memmodelobj.Wp = idmodelobj.M{1}-eye
%   memmodelobj.Wm = idmodelobj.M{2}-eye
%   memmodelobj.w  = idmodelobj.w
%   memmodelobj.fp = idmodelobj.fp


error(CheckSize(idmodelobj,@(x) x.NumPlast==2,'NumPlast==2'));

memmodelobj=SynapseMemoryModel;
memmodelobj.Wp = idmodelobj.M{1}-eye(idmodelobj.NumStates);
memmodelobj.Wm = idmodelobj.M{2}-eye(idmodelobj.NumStates);
memmodelobj.w  = idmodelobj.w;
memmodelobj.fp = idmodelobj.fp;

end

