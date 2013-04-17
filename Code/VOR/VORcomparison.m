function VORcomparison( Wp,Wm,w,t,df,tchange,WpM,WmM,varargin )
%ph=VOREXPT(Wp,Wm,w,t,df,tchange,...) Simulation of VOR expt
%   T = time values
%   WP = potentiation transition rates
%   WM = depression transition rates
%   DF = Change in fraction of potentiation transitions
%   w = Weights of states (+/-1)

Parent=gca;
altfig=[];
fig3=[];
Superpose=false;
axFontSize=40;
FontSize=20;
EqFontSize=20;
txFontSize=20;
Pooled=false;
mutant='D^bK^b-/-';
red=[192 80 77]/256;
varargin=assignApplicable(varargin);

if isempty(altfig)
    altaxWT=[];
    altaxKn=[];
else
    if isscalar(altfig)
        clf(altfig);
        altaxWT=subplot(2,1,1,'Parent',altfig);
        altaxKn=subplot(2,1,2,'Parent',altfig);
        txFontSize=10;
    else
        clf(altfig(1));
        clf(altfig(2));
        altaxWT=axes('Parent',altfig(1));
        altaxKn=axes('Parent',altfig(2));
        embiggen(altaxWT,EqFontSize);
        embiggen(altaxKn,EqFontSize);
    end
end

if Superpose
%     phWT=VORsuperpose(Wp,Wm,w,t,df,'Color','k','Parent',Parent,'altax',altaxWT,varargin{:});
%     phKn=VORsuperpose(WpM,WmM,w,t,df,'Color',red,'Parent',Parent,'altax',altaxKn,varargin{:});
    phWT=VORexptSuperpose(Wp,Wm,w,t,df,tchange,'Color','k','Parent',Parent,varargin{:});
    phKn=VORexptSuperpose(WpM,WmM,w,t,df,tchange,'Color',red,'Parent',Parent,varargin{:});
else
    [phWT,P_WT_nopre,P_WT_pre]=VORexpt(Wp,Wm,w,t,df,tchange,'Color','k','Parent',Parent,varargin{:});
    [phKn,P_Kn_nopre,P_Kn_pre]=VORexpt(WpM,WmM,w,t,df,tchange,'Color',red,'Parent',Parent,varargin{:});
end%if 

PretrainFactor=1;
varargin=assignApplicable(varargin);


if ~isempty(altfig)
    EqProbPlot(Wp,Wm,df,'Parent',altaxWT,'Pooled',Pooled,varargin{:});
    EqProbPlot(WpM,WmM,df,'Parent',altaxKn,'Pooled',Pooled,varargin{:});
%     plot(altax,[EqProb(fp*Wp+(1-fp)*Wm); EqProb((fp+df)*Wp+(1-fp-df)*Wm); EqProb((fp-df)*Wp+(1-fp+df)*Wm)]');
%     xlabel(altax,'State');
%     ylabel(altax,'Equilibrium prob.');
%     legend(altax,{'Untrained','Gain increase','Gain decrease'},'Location','Best');
end

axes(Parent);
embiggen(Parent,FontSize);
xlabel('Training time','FontSize',axFontSize)
ylabel('Learning (-\Delta mean w)','FontSize',axFontSize)

if Superpose
    legend([phWT(1);phKn(1);phWT(2);phKn(2)],{'WT w/o pre';[mutant ' w/o pre'];'WT w/ pre';[mutant ' w/ pre']},'Location','Best')
else
    yl=ylim;
    xl=xlim;
    line(tchange*[1 1],[yl(1) 0],'LineStyle',':','Color','k');
    [x,y]=dsxy2figxy([xl(1) xl(2)], (0.95*yl(2)+0.05*yl(1))*[1 1]);
    annotation('doublearrow',x,y);
    pos=dsxy2figxy([0.6*xl(1)+0.4*xl(2) (0.9*yl(2)+0.1*yl(1)) 0.2*(xl(2)-xl(1)) 0.05*(yl(2)-yl(1))]);
    annotation('textbox',pos,'String','training','VerticalAlignment','top','LineStyle','none','HorizontalAlignment','center','FontSize',txFontSize);
    [x,y]=dsxy2figxy([xl(1) tchange], (0.95*yl(1)+0.05*yl(2))*[1 1]);
    annotation('doublearrow',x,y);
    pos=dsxy2figxy([0.6*xl(1)+0.4*tchange (0.95*yl(1)+0.05*yl(2)) 0.2*(tchange-xl(1)) 0.05*(yl(2)-yl(1))]);
    annotation('textbox',pos,'String','pre-training','VerticalAlignment','bottom','LineStyle','none','HorizontalAlignment','center','FontSize',txFontSize);
    [x,y]=dsxy2figxy([tchange xl(2)], (0.95*yl(1)+0.05*yl(2))*[1 1]);
    annotation('doublearrow',x,y);
    pos=dsxy2figxy([0.6*tchange+0.4*xl(2) (0.95*yl(1)+0.05*yl(2)) 0.2*(xl(2)-tchange) 0.05*(yl(2)-yl(1))]);
    annotation('textbox',pos,'String','training','VerticalAlignment','bottom','LineStyle','none','HorizontalAlignment','center','FontSize',txFontSize);
    legend([phWT(1);phKn(1)],{'WT';'D^bK^b-/-'},'Location','Best')
    drawnow;
    
    if isscalar(fig3)
        clf(fig3);
        h=subplot(2,2,1,'Parent',fig3);
        ProbEvol(P_WT_nopre,t,'WT no pre-training','Parent',h,'Pooled',Pooled);
        h=subplot(2,2,2,'Parent',fig3);
        ProbEvol(P_WT_pre,t,'WT pre-training','Parent',h,'Pooled',Pooled);
        h=subplot(2,2,3,'Parent',fig3);
        ProbEvol(P_Kn_nopre,t,[ mutant ' no pre-training'],'Parent',h,'Pooled',Pooled);
        h=subplot(2,2,4,'Parent',fig3);
        ProbEvol(P_Kn_pre,t,[mutant ' pre-training'],'Parent',h,'Pooled',Pooled);
    elseif numel(fig3)==4
        h=axes('Parent',fig3(1));
        ProbEvol(P_WT_nopre,t,'WT no pre-training','Parent',h,'FontSize',FontSize,'axFontSize',axFontSize,'Pooled',Pooled);
%         embiggen(h,FontSize);
        h=axes('Parent',fig3(2));
        ProbEvol(P_WT_pre,t,'WT pre-training','Parent',h,'FontSize',FontSize,'axFontSize',axFontSize,'Pooled',Pooled);
%         embiggen(h,FontSize);
        h=axes('Parent',fig3(3));
        ProbEvol(P_Kn_nopre,t,[ mutant ' no pre-training'],'Parent',h,'FontSize',FontSize,'axFontSize',axFontSize,'Pooled',Pooled);
%         embiggen(h,FontSize);
        h=axes('Parent',fig3(4));
        ProbEvol(P_Kn_pre,t,[mutant ' pre-training'],'Parent',h,'FontSize',FontSize,'axFontSize',axFontSize,'Pooled',Pooled);
%         embiggen(h,FontSize);
    end
end

if ~isempty(altfig)
    title(altaxWT,'WT');
    title(altaxKn,mutant);
end





end

