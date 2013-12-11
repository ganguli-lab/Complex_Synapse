function [ q ] = GetEnc( obj )
%q=SynapseMemoryModel.GETENC memory encoding matrix
%   q = Wp - Wm
%   WP = potentiation transition rates
%   WM = depression transition rates
%   FP = Fraction of potentiation transitions

q = obj.Wp-obj.Wm;

end

