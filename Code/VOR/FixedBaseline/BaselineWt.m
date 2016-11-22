function [ mw,gradp,gradm ] = BaselineWt( builder_h, n, pot, dep, fpNorm, grad_h )
%dw=BASELINEWT(builder_h,pot_WT,dep_WT,pot_KO,dep_KO) mean baseline synaptic weight.
%   mw = difference in mean baseline synaptic weight
%   builder_h = function handle that builds Wp, Wm and w
%   n      = number of states
%   pot    = parameter for potentiation
%   dep    = parameter for depression
%   fpNorm = fraction of potentiation events at baseline


[Wp,~,w] = builder_h(n, pot);
[~,Wm,~] = builder_h(n, dep);

model = SynapseMemoryModel('Wp',Wp,'Wm',Wm,'fp',fpNorm,'w',w);

mw = model.MeanWeight;

if exist('grad_h','var') && ~isempty(grad_h)
    [Zinv,piv] = model.GetZinv;
    p = piv / Zinv;
    dWp = grad_h(pot, n);
    [~,dWm] = grad_h(dep, n);
    gradp = (p * dWp / Zinv) * model.w * model.fp;
    gradm = (p * dWm / Zinv) * model.w * (1 - model.fp);
end

end

