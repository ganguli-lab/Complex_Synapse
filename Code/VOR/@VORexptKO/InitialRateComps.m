function [ comps ] = InitialRateComps( obj )
%comps=VORexptKO.INITIALRATECOMPS Comparison of initial learning rates in all 4
%cases, no diagonal comparisons
%   obj   = VORexptKO object
%   comps = learning rate differences: [left bottom right top]  
%           +ve = same sign as real experiment
%           -ve = opposite sign to real experiment
%



WT_nopre = obj.nopre.InitialLearningRate(obj.WT);
KO_nopre = obj.nopre.InitialLearningRate(obj.KO);
WT_pre = obj.withpre.InitialLearningRate(obj.WT);
KO_pre = obj.withpre.InitialLearningRate(obj.KO);

comps = [WT_nopre - WT_pre,...
    KO_pre - WT_pre,...
    KO_pre - KO_nopre,...
    WT_nopre - KO_nopre];

end

