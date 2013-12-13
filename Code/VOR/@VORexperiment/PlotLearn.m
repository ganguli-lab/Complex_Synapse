function PlotLearn( obj,varargin )
%VORexperiment.PLOTLEARN plot learning curves during taining
%after any pre-training for WT/KO

Parent=gca;
varargin=assignApplicable(varargin);
cla(Parent);

dt=obj.nopre.tTrain(end)/obj.numpts;

[S,~,t]=obj.nopre.LearningCurve(obj.WT,dt);
ph(1)=plot(t,S(1)-S,'Color',obj.WTcolor,'Parent',Parent,varargin{:});
hold(Parent,'on');
[S]=obj.nopre.LearningCurve(obj.KO,dt);
ph(2)=plot(t,S(1)-S,'Color',obj.KOcolor,'Parent',Parent,varargin{:});
[S]=obj.withpre.LearningCurve(obj.WT,dt);
ph(3)=plot(t,S(1)-S,'Color',obj.WTcolor,'Parent',Parent,varargin{:});
[S]=obj.withpre.LearningCurve(obj.KO,dt);
ph(4)=plot(t,S(1)-S,'Color',obj.KOcolor,'Parent',Parent,varargin{:});

axes(Parent);
embiggen(Parent,obj.FontSize);
xlabel(Parent,'Training time','FontSize',obj.LabFontSize)
ylabel(Parent,'Learning (-\Delta mean w)','FontSize',obj.LabFontSize)

    tchange=obj.withpre.tTrain(1);
    yl=ylim(Parent);
    xl=xlim(Parent);
    line(tchange*[1 1],[yl(1) 0],'LineStyle',':','Color','k','Parent',Parent);
    [x,y]=dsxy2figxy(Parent,[xl(1) xl(2)], (0.95*yl(2)+0.05*yl(1))*[1 1]);
    annotation('doublearrow',x,y);
    pos=dsxy2figxy(Parent,[0.6*xl(1)+0.4*xl(2) (0.9*yl(2)+0.1*yl(1)) 0.2*(xl(2)-xl(1)) 0.05*(yl(2)-yl(1))]);
    annotation('textbox',pos,'String','training','VerticalAlignment','top','LineStyle','none',...
        'HorizontalAlignment','center','FontSize',obj.txFontSize);
    [x,y]=dsxy2figxy(Parent,[xl(1) tchange], (0.95*yl(1)+0.05*yl(2))*[1 1]);
    annotation('doublearrow',x,y);
    pos=dsxy2figxy(Parent,[0.6*xl(1)+0.4*tchange (0.95*yl(1)+0.05*yl(2)) 0.2*(tchange-xl(1)) 0.05*(yl(2)-yl(1))]);
    annotation('textbox',pos,'String','pre-training','VerticalAlignment','bottom','LineStyle','none',...
        'HorizontalAlignment','center','FontSize',obj.txFontSize);
    [x,y]=dsxy2figxy(Parent,[tchange xl(2)], (0.95*yl(1)+0.05*yl(2))*[1 1]);
    annotation('doublearrow',x,y);
    pos=dsxy2figxy(Parent,[0.6*tchange+0.4*xl(2) (0.95*yl(1)+0.05*yl(2)) 0.2*(xl(2)-tchange) 0.05*(yl(2)-yl(1))]);
    annotation('textbox',pos,'String','training','VerticalAlignment','bottom','LineStyle','none',...
        'HorizontalAlignment','center','FontSize',obj.txFontSize);
    legend(Parent,[ph(1);ph(2)],{obj.WTlabel; obj.KOlabel},'Location','Best')
    drawnow;

end

