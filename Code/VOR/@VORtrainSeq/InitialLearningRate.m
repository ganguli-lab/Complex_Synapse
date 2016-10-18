function [ rate ] = InitialLearningRate( obj, modelobj )
%rate=obj.INITIALLEARNINGRATE(modelobj) Initial learning rate
%   Measures learning rate at begining of last epoch
%   obj      = VORtrainSeq object
%   modelobj = SynapseMemoryModel object
%   rate     = initial earning rate: rate of decrease of mean weight
%            = - p * W * w
%   where: p = steady state distribution of W from 2nd last epoch
%          W = transition rate matrix of last epoch
%            = fp * Wp + (1-fp) * Wm 
%          w = synaptic weight vector
%          fp from obj, 
%          Wp,Wm,w from modelobj

modelobj = modelobj.setFp(obj.fps(end-1));
p = modelobj.EqProb();

modelobj = modelobj.setFp(obj.fps(end));
W = modelobj.GetWf();

rate = - p * W * modelobj.w;

end

