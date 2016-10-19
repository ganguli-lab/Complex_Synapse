function EqProbPlots( obj,fh,varargin )
%VORexptKO.EQPROBPLOTS(fh) plot equilibrium distributions for baseline,
%gain-increase, gain-decrease for WT/KO
%   fh = figure handle(s) for plots


n=(0.5:obj.WT.NumStates+0.5)';

if obj.pooled
    n=n-1;
    xlab='Number potentiated';
else
    xlab='State';
end

if isscalar(fh)
    clf(fh);
    h(1)=subplot(1,2,1,'Parent',fh,'FontSize',obj.EqFontSize);
    h(2)=subplot(1,2,2,'Parent',fh,'FontSize',obj.EqFontSize);
elseif numel(fh)==2
    clf(fh(1));clf(fh(2));
    h(1)=axes('Parent',fh(1),'FontSize',obj.EqFontSize);
    h(2)=axes('Parent',fh(2),'FontSize',obj.EqFontSize);
else
    error('Need 1 or 2 figure handles');
end
    
PlotEqProbs(obj.WT,obj.WTlabel,h(1));
PlotEqProbs(obj.KO,obj.KOlabel,h(2));




    function PlotEqProbs(modelobj,titletext,Parent)
        modelobj=modelobj.setFp(obj.withpre.fps(1));
        p=modelobj.EqProb;
        modelobj=modelobj.setFp(obj.withpre.fps(3));
        p=[p;modelobj.EqProb];
        modelobj=modelobj.setFp(obj.withpre.fps(2));
        p=[p;modelobj.EqProb]';
        %if we're using stairs:
        p=[p;p(end,:)];

        stairs(Parent,n,p,varargin{:});
        % plot(Parent,n(1:end-1)+0.5,p,varargin{:});
        % bar(Parent,n(1:end-1)+0.5,p,varargin{:});
        xlim(Parent,[n(1) n(end)]);
        set(Parent,'XTick',n(1:end-1)+0.5);
        xlabel(Parent,xlab);
        ylabel(Parent,'Equilibrium probability');
        title(Parent,titletext);
        legend(Parent,{'Untrained','Gain increase','Gain decrease'},'Location','Best');
    end


end

