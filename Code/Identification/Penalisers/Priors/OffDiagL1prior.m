function [ logprior ] = OffDiagL1prior( modelobj,options )
%logprior=OFFDIAGL1PRIOR(modelobj,options) logarithm of unnormalised
%off-diagonal L^1 prior evaluated on SynapeIdModel modelobj
%   logprior = logarithm of unnormalised prior
%   modelobj = SynapeIdModel
%   options  = struct of options (see SynapseOptimset)


logprior=0;
for i=1:modelobj.NumPlast
    logprior = logprior - sum(sum(modelobj.M{i})) - trace(modelobj.M{i});
end
logprior = logprior * options.Penalty;

end

