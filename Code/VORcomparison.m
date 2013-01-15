function VORcomparison( Wp,Wm,w,t,df,tchange,WpM,WmM,varargin )
%ph=VOREXPT(Wp,Wm,w,t,df,tchange,...) Simulation of VOR expt
%   T = time values
%   WP = potentiation transition rates
%   WM = depression transition rates
%   DF = Change in fraction of potentiation transitions
%   w = Weights of states (+/-1)

phWT=VORexpt(Wp,Wm,w,t,df,tchange,'LinSpec','k',varargin{:});
phKn=VORexpt(WpM,WmM,w,t,df,tchange,'LinSpec','r',varargin{:});
xlabel('Training time')
ylabel('Learning (change in mean synaptic weight)')
yl=ylim;
xl=xlim;
[x,y]=dsxy2figxy([xl(1) xl(2)], 0.95*yl(2)*[1 1]);
annotation('doublearrow',x,y);
pos=dsxy2figxy([0.6*xl(1)+0.4*xl(2) 0.9*yl(2) 0.2*(xl(2)-xl(1)) 0.05*yl(2)]);
annotation('textbox',pos,'String','training','VerticalAlignment','top','LineStyle','none','HorizontalAlignment','center');
[x,y]=dsxy2figxy([xl(1) tchange], 0.95*yl(1)*[1 1]);
annotation('doublearrow',x,y);
pos=dsxy2figxy([0.6*xl(1)+0.4*tchange 0.95*yl(1) 0.2*(tchange-xl(1)) abs(0.05*yl(1))]);
annotation('textbox',pos,'String','pre-training','VerticalAlignment','bottom','LineStyle','none','HorizontalAlignment','center');
[x,y]=dsxy2figxy([tchange xl(2)], 0.95*yl(1)*[1 1]);
annotation('doublearrow',x,y);
pos=dsxy2figxy([0.6*tchange+0.4*xl(2) 0.95*yl(1) 0.2*(xl(2)-tchange) abs(0.05*yl(1))]);
annotation('textbox',pos,'String','training','VerticalAlignment','bottom','LineStyle','none','HorizontalAlignment','center');
legend([phWT(1);phKn(1)],{'WT';'D^bK^b-/-'},'Location','Best')
embiggen

end

