function [ comp ] = InitRComp_top( obj )
%comps=obj.INITRCOMP_top Comparison of initial learning rates in top 2
%cases
%   obj  = VORexperiment object
%   comp = learning rate differences: [left bottom right top]  
%          +ve = same sign as real experiment
%          -ve = opposite sign to real experiment
%



WT_nopre = obj.nopre.InitialLearningRate(obj.WT);
KO_nopre = obj.nopre.InitialLearningRate(obj.KO);

comp = WT_nopre - KO_nopre;

end

