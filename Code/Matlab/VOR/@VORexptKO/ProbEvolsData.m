function [ P_WT_nopre,P_KO_nopre,P_WT_pre,P_KO_pre,t ] = ProbEvolsData( obj )
%[P_WT_nopre,P_KO_nopre,P_WT_pre,P_KO_pre,t]=VORexptKO.PROBEVOLSDATA 
%evolution of distributions during taining and pre-training for WT/KO
%   P_XX_yyy = distribution across synaptic states (2nd ind) as a function
%               of time (1st ind)
%   XX       = WT/KO
%   yyy      = nopre/pre
%   t        = time


dt=obj.nopre.tTrain(end)/obj.numpts;

[~,P_WT_nopre,t]=obj.nopre.LearningCurve(obj.WT,dt);
[~,P_KO_nopre]=obj.nopre.LearningCurve(obj.KO,dt);
[~,P_WT_pre]=obj.withpre.LearningCurve(obj.WT,dt);
[~,P_KO_pre]=obj.withpre.LearningCurve(obj.WT,dt);


end

