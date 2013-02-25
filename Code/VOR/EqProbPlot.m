function EqProbPlot( Wp,Wm,df,varargin )
%EQPROBPLOT(Wp,Wm,df,...) Summary of this function goes here
%   Detailed explanation goes here

fp=0.5;
Parent=gca;
Pooled=false;

n=(0.5:length(Wp)+0.5)';
p=[EqProb(fp*Wp+(1-fp)*Wm); EqProb((fp+df)*Wp+(1-fp-df)*Wm); EqProb((fp-df)*Wp+(1-fp+df)*Wm)]';
p=[p;p(end,:)];

varargin=assignApplicable(varargin);

if Pooled
    n=n-1;
    xlab='Number potentiated';
else
    xlab='State';
end

stairs(Parent,n,p,varargin{:});
% plot(Parent,[EqProb(fp*Wp+(1-fp)*Wm); EqProb((fp+df)*Wp+(1-fp-df)*Wm); EqProb((fp-df)*Wp+(1-fp+df)*Wm)]',varargin{:});
xlabel(Parent,xlab);
xlim(Parent,[n(1) n(end)]);
set(Parent,'XTick',n(1:end-1)+0.5);
ylabel(Parent,'Equilibrium probability');
legend(Parent,{'Untrained','Gain increase','Gain decrease'},'Location','Best');

end

