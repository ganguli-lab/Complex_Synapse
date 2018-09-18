function ProbEvols( obj,fh,varargin )
%VORexptKO.PROBEVOLS(fh) plot evolution of distributions during taining
%and pre-training for WT/KO
%   fh = figure handle(s) for plots


n=1:obj.WT.NumStates;
if obj.pooled
    n=n-1;
    xlab='Number potentiated';
else
    xlab='State';
end

[P_WT_nopre,P_KO_nopre,P_WT_pre,P_KO_pre,t]=obj.ProbEvolsData;



if isscalar(fh)
    clf(fh);
    h(1)=subplot(2,2,1,'Parent',fh);
    h(2)=subplot(2,2,2,'Parent',fh);
    h(3)=subplot(2,2,3,'Parent',fh);
    h(4)=subplot(2,2,4,'Parent',fh);
elseif numel(fh)==4
    clf(fh(1));clf(fh(2));clf(fh(3));clf(fh(4));
    h(1)=axes('Parent',fh(1));
    h(2)=axes('Parent',fh(2));
    h(3)=axes('Parent',fh(3));
    h(4)=axes('Parent',fh(4));
else
    error('Need 1 or 4 figure handles');
end
    
PlotProbEvol(P_WT_nopre,[obj.WTlabel ' ' obj.noprelabel],h(1));
PlotProbEvol(P_KO_nopre,[obj.KOlabel ' ' obj.noprelabel],h(2));
PlotProbEvol(P_WT_pre,[obj.WTlabel ' ' obj.withprelabel],h(3));
PlotProbEvol(P_KO_pre,[obj.KOlabel ' ' obj.withprelabel],h(4));

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

