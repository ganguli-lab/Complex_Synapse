function [ Wf ] = GetWf( obj )
%Wf=SynapseMemoryModel.GETWF forgetting matrix
%   Wf = fp*Wp+fm*Wm
%   WP = potentiation transition rates
%   WM = depression transition rates
%   FP = Fraction of potentiation transitions

Wf = obj.fp*obj.Wp+(1-obj.fp)*obj.Wm;

end

