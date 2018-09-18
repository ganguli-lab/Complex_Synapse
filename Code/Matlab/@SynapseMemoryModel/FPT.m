function [ T ] = FPT( obj )
%T=SynapseMemoryModel.FPT Calculate off diagonal mean first passage times
%   W = transition rates


E=ones(size(obj.Wp));
p=obj.EqProb;
Z=inv(obj.GetZinv);

T=(E*diag(diag(Z)) - Z)/diag(p);

end

