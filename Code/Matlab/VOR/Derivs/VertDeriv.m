function [ vder ] = VertDeriv( modelobj, fpInc )
%hder=HORIZDERIV(modelobj,fpInc) Derivative of initial learning rate wrt
%fpNorm - baseline fraction of potentiation events
%   vder = gradient, can have both signs for a viable model
%   modelobj = SynapseMemoryModel object
%   fpInc = fraction of potentiation events, fp, during training


[Zinv,piv] = modelobj.GetZinv;
ZQ = Zinv \ modelobj.GetEnc;
% eqp = piv / Zinv;

dPhidfp = (modelobj.fp - fpInc) * ZQ * ZQ  + ZQ;
vder = piv * dPhidfp * modelobj.w;



end

