function [ logprior ] = Lhalfprior( modelobj,options )
%logprior=LHALFPRIOR(modelobj,options) logarithm of unnormalised
%off-diagonal L^1 prior evaluated on SynapeIdModel modelobj 
%   logprior = logarithm of unnormalised prior
%   modelobj = SynapeIdModel
%   options  = struct of options (see SynapseOptimset)


logprior=0;
for i=1:modelobj.NumPlast
    logprior = logprior - sum(modelobj.M{i}(:).^(1/2));
end
logprior = logprior * options.Penalty;

end

