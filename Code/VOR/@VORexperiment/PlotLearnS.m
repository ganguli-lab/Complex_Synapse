function PlotLearnS( obj,varargin )
%VORexperiment.PLOTLEARNS plot learning curves during taining
%and pre-training for WT/KO


persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='VORexperiment.PlotLearnS';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addOptional('extraArgs',{},@(x)validateattributes(x,{'cell'},{},'VORexperiment.PlotLearnS','extraArgs',2));
    p.addParameter('Parent',gca,@(x)validateattributes(x,{'numeric'},{'scalar'},'VORexperiment.PlotLearnS','Mlabels'));
%     p.addParameter('Normalise',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
end
p.parse(varargin{:});
r=p.Results;

cla(r.Parent);

tchange=obj.withpre.tTrain(1);

dt=diff(obj.withpre.tTrain)/obj.numpts;


[S,~,t]=obj.nopre.LearningCurve(obj.WT,dt);
T=t(end)-tchange;
ph(1)=plot(t(t<=T),S(1)-S(t<=T),'Color',obj.WTcolor,'LineStyle',obj.noprestyle,'Parent',r.Parent,r.extraArgs{:});
hold(r.Parent,'on');
[S]=obj.withpre.LearningCurve(obj.WT,dt);
S(t<tchange)=[];
ph(3)=plot(t(t>=tchange)-tchange,S(1)-S,'Color',obj.WTcolor,'LineStyle',obj.withprestyle,'Parent',r.Parent,r.extraArgs{:});

[S,~,t]=obj.nopre.LearningCurve(obj.KO,dt);
T=t(end)-tchange;
ph(2)=plot(t(t<=T),S(1)-S(t<=T),'Color',obj.KOcolor,'LineStyle',obj.noprestyle,'Parent',r.Parent,r.extraArgs{:});
[S]=obj.withpre.LearningCurve(obj.KO,dt);
S(t<tchange)=[];
ph(4)=plot(t(t>=tchange)-tchange,S(1)-S,'Color',obj.KOcolor,'LineStyle',obj.withprestyle,'Parent',r.Parent,r.extraArgs{:});



set(r.Parent,'FontSize',obj.FontSize);
xlabel(r.Parent,'Training time','FontSize',obj.LabFontSize)
ylabel(r.Parent,'Learning (-\Delta mean w)','FontSize',obj.LabFontSize)
legend(r.Parent,ph,{[obj.WTlabel ' ' obj.noprelabel];[obj.KOlabel ' ' obj.noprelabel];...
    [obj.WTlabel ' ' obj.withprelabel];[obj.KOlabel ' ' obj.withprelabel]},'Location','Best')






end

