function PlotLearnS( obj,varargin )
%VORexperiment.PLOTLEARNS plot learning curves during taining
%after any pre-training for WT/KO


persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='VORexperiment.PlotLearnS';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addParameter('Parent',gca,@(x)validateattributes(x,{'numeric','matlab.graphics.axis.Axes'},{'scalar'},'VORexperiment.PlotLearnS','Parent'));
end
p.parse(varargin{:});
r=p.Results;

cla(r.Parent);


St=obj.LearnSdata();

ph(1)=plot(St.t,St.WTnopre,'Color',obj.noprecolor,'LineStyle',obj.noprestyle,'Parent',r.Parent,p.Unmatched);
hold(r.Parent,'on');
ph(2)=plot(St.t,St.WTwithpre,'Color',obj.withprecolor,'LineStyle',obj.withprestyle,'Parent',r.Parent,p.Unmatched);




set(r.Parent,'FontSize',obj.FontSize);
xlabel(r.Parent,'Training time','FontSize',obj.LabFontSize)
ylabel(r.Parent,'Learning (-\Delta mean w)','FontSize',obj.LabFontSize)
legend(r.Parent,ph,{obj.noprelabel; obj.withprelabel},'Location','NorthWest','FontSize',obj.LegFontSize)






end

