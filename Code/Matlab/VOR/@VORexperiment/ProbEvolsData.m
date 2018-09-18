function [ P_WT_nopre,P_WT_pre,t ] = ProbEvolsData( obj )
%[P_WT_nopre,P_WT_pre,t]=VORexperiment.PROBEVOLSDATA 
%evolution of distributions during taining and pre-training for WT/KO
%   P_XX_yyy = distribution across synaptic states (2nd ind) as a function
%               of time (1st ind)
%   XX       = WT
%   yyy      = nopre/pre
%   t        = time


dt=obj.nopre.tTrain(end)/obj.numpts;

[~,P_WT_nopre,t]=obj.nopre.LearningCurve(obj.WT,dt);
[~,P_WT_pre]=obj.withpre.LearningCurve(obj.WT,dt);


end

