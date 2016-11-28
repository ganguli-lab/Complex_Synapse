function PlotLearnS( obj,varargin )
%VORexptKO.PLOTLEARNS plot learning curves during taining
%after any pre-training for WT/KO


persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='VORexptKO.PlotLearnS';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addParameter('Parent',gca,@(x)validateattributes(x,{'numeric','matlab.graphics.axis.Axes'},{'scalar'},'VORexptKO.PlotLearnS','Parent'));
    p.addParameter('all',true,@(x)validateattributes(x,{'numeric','logical'},{'scalar'},'VORexptKO.PlotLearnS','all'));
end
p.parse(varargin{:});
r=p.Results;

if ~isgraphics(r.Parent)
    r.Parent=gca;
end

cla(r.Parent);


St=obj.LearnSdata();

if r.all
    ph(1)=plot(St.t,St.WTnopre,'Color',obj.WTcolor,'LineStyle',obj.noprestyle,'Parent',r.Parent,p.Unmatched);
    hold(r.Parent,'on');
    ph(3)=plot(St.t,St.WTwithpre,'Color',obj.WTcolor,'LineStyle',obj.withprestyle,'Parent',r.Parent,p.Unmatched);

    ph(2)=plot(St.t,St.KOnopre,'Color',obj.KOcolor,'LineStyle',obj.noprestyle,'Parent',r.Parent,p.Unmatched);
    ph(4)=plot(St.t,St.KOwithpre,'Color',obj.KOcolor,'LineStyle',obj.withprestyle,'Parent',r.Parent,p.Unmatched);
    
    lab = {[obj.WTlabel ' ' obj.noprelabel];[obj.KOlabel ' ' obj.noprelabel];...
    [obj.WTlabel ' ' obj.withprelabel];[obj.KOlabel ' ' obj.withprelabel]};
else
    ph(1)=plot(St.t,St.WTwithpre,'Color',obj.WTcolor,'LineStyle',obj.withprestyle,'Parent',r.Parent,p.Unmatched);
    hold(r.Parent,'on');
    ph(2)=plot(St.t,St.KOwithpre,'Color',obj.KOcolor,'LineStyle',obj.withprestyle,'Parent',r.Parent,p.Unmatched);

    lab = {[obj.WTlabel ' ' obj.withprelabel];[obj.KOlabel ' ' obj.withprelabel]};
end

set(r.Parent,'FontSize',obj.FontSize);
xlabel(r.Parent,'Training time','FontSize',obj.LabFontSize)
ylabel(r.Parent,'Learning (-\Delta mean w)','FontSize',obj.LabFontSize)
legend(r.Parent,ph,lab,'Location','Best')






end

