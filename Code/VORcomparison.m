function VORcomparison( Wp,Wm,w,t,df,tchange,WpM,WmM,varargin )
%ph=VOREXPT(Wp,Wm,w,t,df,tchange,...) Simulation of VOR expt
%   T = time values
%   WP = potentiation transition rates
%   WM = depression transition rates
%   DF = Change in fraction of potentiation transitions
%   w = Weights of states (+/-1)

phWT=VORexpt(Wp,Wm,w,t,df,tchange,varargin{:});
phKn=VORexpt(WpM,WmM,w,t,df,tchange,'LinSpec','g',varargin{:});
legend([phWT(1);phKn(1)],{'WT';'Knockout'},'Location','Best')
xlabel('Training time')
ylabel('Learning (mean synaptic weight)')
embiggen

end

