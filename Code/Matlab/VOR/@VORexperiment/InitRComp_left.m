function [ comp ] = InitRComp_left( obj )
%comps=obj.INITRCOMP_left Comparison of initial learning rates in left 2
%cases
%   obj  = VORexperiment object
%   comp = learning rate differences: [left bottom right top]  
%          +ve = same sign as real experiment
%          -ve = opposite sign to real experiment
%



WT_nopre = obj.nopre.InitialLearningRate(obj.WT);
WT_pre = obj.withpre.InitialLearningRate(obj.WT);

comp = WT_nopre - WT_pre;

end

