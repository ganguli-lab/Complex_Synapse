function [ Zinv,piv ] = GetZinv( obj,piv )
%[Zinv,piv]=SynapseMemoryModel.GETZINV inverse of fundamental matrix
%   Zinv = ev*piv - Wf
%   ev = column vector of ones
%   piv = arbitrary row vec, s.t. piv*ev~=0 (default ev')

ev=ones(size(obj.w));

existsAndDefault('piv',ev');

Zinv = ev*piv - obj.GetWf;


end

