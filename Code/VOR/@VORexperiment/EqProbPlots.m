function EqProbPlots( obj,fh,varargin )
%VORexperiment.EQPROBPLOTS(fh) plot equilibrium distributions for baseline,
%gain-increase, gain-decrease for WT
%   fh = figure handle(s) for plots


n=(0.5:obj.WT.NumStates+0.5)';

if obj.pooled
    n=n-1;
    xlab='Number potentiated';
else
    xlab='State';
end

h=axes('Parent',fh,'FontSize',obj.EqFontSize);

PlotEqProbs(obj.WT,obj.WTlabel,h);




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

