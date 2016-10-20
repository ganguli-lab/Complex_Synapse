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
        leg = {'Untrained'};
        modelobj = modelobj.setFp(obj.withpre.fps(1));
        p = modelobj.EqProb;
        if obj.withpre.numTrain >= 2
            modelobj = modelobj.setFp(obj.withpre.fps(2));
            p = [p; modelobj.EqProb];
            leg{end+1} = 'After pretraining';
        end
        if obj.withpre.numTrain >= 1
            modelobj = modelobj.setFp(obj.withpre.fps(end));
            p = [p; modelobj.EqProb];
            leg{end+1} = 'After training';
        end
        %if we're using stairs:
        p = [p p(:,end)]';

        stairs(Parent, n, p, varargin{:});
        % plot(Parent,n(1:end-1)+0.5,p,varargin{:});
        % bar(Parent,n(1:end-1)+0.5,p,varargin{:});
        xlim(Parent, [n(1) n(end)]);
        set(Parent, 'XTick', n(1:end-1)+0.5);
        xlabel(Parent, xlab);
        ylabel(Parent, 'Equilibrium probability');
        title(Parent, titletext);
        legend(Parent, leg, 'Location', 'Best');
    end


end

