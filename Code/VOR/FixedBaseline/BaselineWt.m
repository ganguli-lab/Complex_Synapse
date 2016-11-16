function [ mw ] = BaselineWt( builder_h, n, pot, dep, fpNorm )
%dw=BASELINEWT(builder_h,pot_WT,dep_WT,pot_KO,dep_KO) mean baseline synaptic weight.
%   mw = difference in mean baseline synaptic weight
%   builder_h = function handle that builds Wp, Wm and w
%   n      = number of states
%   pot = parameter for potentiation
%   dep = parameter for depression
%   fpNorm = fraction of potentiation events at baseline


[Wp,~,w] = builder_h(n, pot);
[~,Wm,~] = builder_h(n, dep);

model = SynapseMemoryModel('Wp',Wp,'Wm',Wm,'fp',fpNorm,'w',w);

mw = model.MeanWeight;

end

