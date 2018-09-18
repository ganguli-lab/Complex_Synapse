function [ dw ] = BaseDiff( builder_h, n, pot_WT, dep_WT, pot_KO, dep_KO, fpNorm )
%dw=BASEDIFF(builder_h,pot_WT,dep_WT,pot_KO,dep_KO) difference in mean
%baseline synaptic weight.
%   dw = difference in mean baseline synaptic weight
%   builder_h = function handle that builds Wp, Wm and w
%   n      = number of states
%   pot_?? = parameter for potentiation
%   dep_?? = parameter for depression
%   ???_WT = wild-type parameter
%   ???_KO = wild-type parameter
%   fpNorm = fraction of potentiation events at baseline


[Wp,~,w] = builder_h(n, pot_WT);
[~,Wm,~] = builder_h(n, dep_WT);
wtmodel = SynapseMemoryModel('Wp',Wp,'Wm',Wm,'fp',fpNorm,'w',w);

[Wp,~,w] = builder_h(n, pot_KO);
[~,Wm,~] = builder_h(n, dep_KO);
komodel = SynapseMemoryModel('Wp',Wp,'Wm',Wm,'fp',fpNorm,'w',w);

dw = wtmodel.MeanWeight - komodel.MeanWeight;

end

