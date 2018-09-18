function [ hder ] = HorizDeriv( modelobj, dWp, dWm, fixWt )
%hder=HORIZDERIV(modelobj,dWp,dWm,fixWt) Derivative of initial learning
%rate wrt depression parameter, divided by Delta f
%   hder = gradient, can have both signs for a viable model
%   modelobj = SynapseMemoryModel object
%   dWP = gradient of potentiation transition rates
%   dWM = gradient of depression transition rates
%   fixWt = boolean, adjust Wp to keep p.w fixed?

if iscell(dWp)
    hder = cell(size(dWp));
    for i = 1:numel(dWp)
        hder{i} = HorizDeriv( modelobj, dWp{i}, dWm{i}, fixWt );
    end%for i

else%if iscell

    [Zinv,piv] = modelobj.GetZinv;
    dQ = modelobj.GetEnc;
    eqp = piv / Zinv;

    dPhidxm = (1-modelobj.fp) * (dWm / Zinv) * dQ  -  dWm;
    hder = eqp * dPhidxm * modelobj.w;

    if fixWt
        dPhidxp = modelobj.fp * (dWp / Zinv) * dQ  + dWp; 
        dwdxp = modelobj.fp * (eqp * dWp / Zinv) * modelobj.w;
        dwdxm = (1-modelobj.fp) * (eqp * dWm / Zinv) * modelobj.w;

        dxpdxm = - dwdxm / dwdxp;
        hder = hder + dxpdxm *  eqp * dPhidxp * modelobj.w;
    end%if fixWt
    
end%if iscell


end

