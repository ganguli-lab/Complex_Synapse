function [ rate ] = InitialLearningRate( obj, modelobj )
%rate=VORtrainSeqDiff.INITIALLEARNINGRATE(modelobj) Initial learning rate
%   Measures learning rate at begining of last epoch
%   obj      = VORtrainSeqDiff object
%   modelobj = SynapseMemoryModel object
%   rate     = initial earning rate: rate of decrease of mean weight
%            = - p * W * w
%   where: p = steady state distribution of W from 2nd last epoch
%          W = transition rate matrix of last epoch
%            = fp * Wp + (1-fp) * Wm 
%          w = synaptic weight vector
%          fp from obj, 
%          Wp,Wm,w from modelobj


rate_rel = obj.VORrel.InitialLearningRate(modelobj);

rate_other = obj.VORcomp.InitialLearningRate(modelobj);

rate = (1 - obj.frac_other) * rate_rel - obj.frac_other * rate_other;


end

