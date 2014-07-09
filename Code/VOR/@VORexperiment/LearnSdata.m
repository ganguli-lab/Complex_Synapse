function St=LearnSdata( obj,varargin )
%VORexperiment.PLOTLEARNS plot learning curves during taining
%and pre-training for WT/KO


tchange=obj.withpre.tTrain(1);

dt=diff(obj.withpre.tTrain)/obj.numpts;


[S,~,t]=obj.nopre.LearningCurve(obj.WT,dt);
T=t(end)-tchange;
% ph(1)=plot(t(t<=T),S(1)-S(t<=T),'Color',obj.WTcolor,'LineStyle',obj.noprestyle,'Parent',Parent,varargin{:});
St.t=t(t<=T);
St.WTnopre=S(1)-S(t<=T);
% hold(Parent,'on');
[S]=obj.withpre.LearningCurve(obj.WT,dt);
S(t<tchange)=[];
% ph(3)=plot(t(t>=tchange)-tchange,S(1)-S,'Color',obj.WTcolor,'LineStyle',obj.withprestyle,'Parent',Parent,varargin{:});
St.WTwithpre=S(1)-S;

[S,~,t]=obj.nopre.LearningCurve(obj.KO,dt);
T=t(end)-tchange;
% ph(2)=plot(t(t<=T),S(1)-S(t<=T),'Color',obj.KOcolor,'LineStyle',obj.noprestyle,'Parent',Parent,varargin{:});
St.KOnopre=S(1)-S(t<=T);
[S]=obj.withpre.LearningCurve(obj.KO,dt);
S(t<tchange)=[];
% ph(4)=plot(t(t>=tchange)-tchange,S(1)-S,'Color',obj.KOcolor,'LineStyle',obj.withprestyle,'Parent',Parent,varargin{:});
St.KOwithpre=S(1)-S;



% embiggen(Parent,obj.FontSize);
% xlabel(Parent,'Training time','FontSize',obj.LabFontSize)
% ylabel(Parent,'Learning (-\Delta mean w)','FontSize',obj.LabFontSize)
% legend(Parent,ph,{[obj.WTlabel ' ' obj.noprelabel];[obj.KOlabel ' ' obj.noprelabel];...
%     [obj.WTlabel ' ' obj.withprelabel];[obj.KOlabel ' ' obj.withprelabel]},'Location','Best')






end

