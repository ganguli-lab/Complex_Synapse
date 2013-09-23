function [ loglike ] = ChunkedHMMloglike( chunks,modelobj,simobj )
%loglike=CHUNKEDHMMLOGLIKE(chunks,modelobj,simobj) likelihood of outputs for
%Hidden-Markov-Model
%   loglike  = log likelihood
%   chunks   = 2-by-K matrix of starts and ends of each chunk.
%   modelobj = SynapseIdModel
%   simobj   = SynapsePlastSeq


loglike=0;
for i=1:size(chunks,2)
    loglike=loglike+HMMloglike(modelobj,simobj.GetRange(chunks(:,i)));
end


end

