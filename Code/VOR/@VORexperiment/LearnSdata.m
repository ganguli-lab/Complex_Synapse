function St=LearnSdata( obj,varargin )
%VORexperiment.PLOTLEARNS plot learning curves during taining
%and pre-training for WT

if obj.withpre.numTrain > 1
%     tchange=obj.withpre.tTrain(1);
    dt=sum(diff(obj.withpre.tTrain))/obj.numpts;
else
%     tchange=0;
    dt=obj.withpre.tTrain/obj.numpts;
end

[S,~,t]=obj.nopre.LearningCurveEnd(obj.WT,dt);
% T=t(end)-tchange;
% ph(1)=plot(t(t<=T),S(1)-S(t<=T),'Color',obj.WTcolor,'LineStyle',obj.noprestyle,'Parent',Parent,varargin{:});
St.t=t;%(t<=T);
St.WTnopre=S(1)-S;%(t<=T);
% hold(Parent,'on');
[S]=obj.withpre.LearningCurveEnd(obj.WT,dt);
% S(t<tchange)=[];
% ph(3)=plot(t(t>=tchange)-tchange,S(1)-S,'Color',obj.WTcolor,'LineStyle',obj.withprestyle,'Parent',Parent,varargin{:});
St.WTwithpre=S(1)-S;









end

