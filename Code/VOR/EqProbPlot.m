function EqProbPlot( Wp,Wm,df,varargin )
%EQPROBPLOT(Wp,Wm,df,...) Summary of this function goes here
%   Detailed explanation goes here

fp=0.5;
Parent=gca;

varargin=assignApplicable(varargin);

plot(Parent,[EqProb(fp*Wp+(1-fp)*Wm); EqProb((fp+df)*Wp+(1-fp-df)*Wm); EqProb((fp-df)*Wp+(1-fp+df)*Wm)]',varargin{:});
xlabel(Parent,'State');
xlim(Parent,[1 length(Wp)]);
ylabel(Parent,'Equilibrium prob.');
legend(Parent,{'Untrained','Gain increase','Gain decrease'},'Location','Best');

end

