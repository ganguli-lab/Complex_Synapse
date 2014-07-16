function PlotLearnS( obj,varargin )
%VORexperiment.PLOTLEARNS plot learning curves during taining
%and pre-training for WT/KO


persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='VORexperiment.PlotLearnS';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addParameter('Parent',gca,@(x)validateattributes(x,{'numeric'},{'scalar'},'VORexperiment.PlotLearnS','Mlabels'));
end
p.parse(varargin{:});
r=p.Results;

cla(r.Parent);

tchange=obj.withpre.tTrain(1);

dt=diff(obj.withpre.tTrain)/obj.numpts;


[S,~,t]=obj.nopre.LearningCurve(obj.WT,dt);
T=t(end)-tchange;
ph(1)=plot(t(t<=T),S(1)-S(t<=T),'Color',obj.WTcolor,'LineStyle',obj.noprestyle,'Parent',r.Parent,p.Unmatched);
hold(r.Parent,'on');
[S]=obj.withpre.LearningCurve(obj.WT,dt);
S(t<tchange)=[];
ph(3)=plot(t(t>=tchange)-tchange,S(1)-S,'Color',obj.WTcolor,'LineStyle',obj.withprestyle,'Parent',r.Parent,p.Unmatched);

[S,~,t]=obj.nopre.LearningCurve(obj.KO,dt);
T=t(end)-tchange;
ph(2)=plot(t(t<=T),S(1)-S(t<=T),'Color',obj.KOcolor,'LineStyle',obj.noprestyle,'Parent',r.Parent,p.Unmatched);
[S]=obj.withpre.LearningCurve(obj.KO,dt);
S(t<tchange)=[];
ph(4)=plot(t(t>=tchange)-tchange,S(1)-S,'Color',obj.KOcolor,'LineStyle',obj.withprestyle,'Parent',r.Parent,p.Unmatched);



set(r.Parent,'FontSize',obj.FontSize);
xlabel(r.Parent,'Training time','FontSize',obj.LabFontSize)
ylabel(r.Parent,'Learning (-\Delta mean w)','FontSize',obj.LabFontSize)
legend(r.Parent,ph,{[obj.WTlabel ' ' obj.noprelabel];[obj.KOlabel ' ' obj.noprelabel];...
    [obj.WTlabel ' ' obj.withprelabel];[obj.KOlabel ' ' obj.withprelabel]},'Location','Best')






end

