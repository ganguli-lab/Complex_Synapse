function [ logprior ] = SynapsePrior( modelobj,options )
%logprior=SYNAPSEPRIOR(modelobj,options) logarithm of unnormalised prior
%evaluated on SynapeIdModel modelobj  
%   logprior = logarithm of unnormalised prior
%   modelobj = SynapeIdModel
%   options  = struct of options (see SynapseOptimset)

priorfn=str2func([options.Penaliser 'prior']);
logprior=priorfn(modelobj,options);

end

