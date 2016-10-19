function ProbEvols( obj,fh,varargin )
%VORexperiment.PROBEVOLS(fh) plot evolution of distributions during taining
%and pre-training for WT/KO
%   fh = figure handle(s) for plots


n=1:obj.WT.NumStates;
if obj.pooled
    n=n-1;
    xlab='Number potentiated';
else
    xlab='State';
end

[P_WT_nopre,P_WT_pre,t]=obj.ProbEvolsData;



if isscalar(fh)
    clf(fh);
    h(1)=subplot(1,2,1,'Parent',fh);
    h(2)=subplot(1,2,2,'Parent',fh);
elseif numel(fh)==2
    clf(fh(1));clf(fh(2))
    h(1)=axes('Parent',fh(1));
    h(2)=axes('Parent',fh(2));
else
    error('Need 1 or 2 figure handles');
end
    
PlotProbEvol(P_WT_nopre,[obj.noprelabel],h(1));
PlotProbEvol(P_WT_pre,[obj.withprelabel],h(2));

drawnow;

    function PlotProbEvol( Pt,titletext,Parent )

    cla(Parent);

    imagesc(t,n,Pt','Parent',Parent,varargin{:});
    ylabel(Parent,xlab,'FontSize',obj.ProbFontSize);
    xlabel(Parent,'Training time','FontSize',obj.ProbFontSize);
    title(Parent,titletext,'FontSize',obj.ProbFontSize);
    ch=colorbar('peer',Parent);
    % set(get(h,'YLabel'),'String','Probability','Rotation',270,'VerticalAlignment','bottom','FontSize',FontSize);
    colorbarlabel(ch,'Probability','FontSize',obj.ProbFontSize);

    end


end

