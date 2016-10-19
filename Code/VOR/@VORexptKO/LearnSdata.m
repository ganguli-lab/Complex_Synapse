function St=LearnSdata( obj,varargin )
%VORexptKO.PLOTLEARNS plot learning curves during taining
%and pre-training for WT/KO


tchange=obj.withpre.tTrain(1);

dt=diff(obj.withpre.tTrain)/obj.numpts;


[S,~,t]=obj.nopre.LearningCurve(obj.WT,dt);
T=t(end)-tchange;
St.t=t(t<=T);
St.WTnopre=S(1)-S(t<=T);

[S]=obj.withpre.LearningCurve(obj.WT,dt);
S(t<tchange)=[];
St.WTwithpre=S(1)-S;

[S,~,t]=obj.nopre.LearningCurve(obj.KO,dt);
T=t(end)-tchange;
St.KOnopre=S(1)-S(t<=T);

[S]=obj.withpre.LearningCurve(obj.KO,dt);
S(t<tchange)=[];
St.KOwithpre=S(1)-S;









end

