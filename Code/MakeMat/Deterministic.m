function [ Wp,Wm,w ] = Deterministic( targets )
%[Wp,Wm,w]=DETERMINISTIC(targets) deterministic synapse
%   Potentiation goes from i to targets(i) with probability 1.
%   Depression is flip of potentiation (i -> n-i+1)
%   WP = potentiation transition rates
%   WM = depression transition rates
%   w = Weights of states [-1;...;-1;1;...;1].

error(CheckSize(targets,@isrow));

Wp=zeros(length(targets));
inds=sub2ind(size(Wp),1:length(targets),targets);
Wp(inds)=1;
Wp=StochastifyC(Wp);

Wm=rot90(Wp,2);

if nargout>2
    w=BinaryWeights(length(targets));
end

end

