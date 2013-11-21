function OptL_GUI( fp,w,Aenv,srange )
%OPTL_GUI Summary of this function goes here
%   Detailed explanation goes here


%-------------------------------------------------------------------------
%Create figure
% figure1 = figure('Units','normalized','OuterPosition',[0 0.03 1 0.97]);
figure1 = figure('Units','normalized','WindowStyle','docked');

%Create panels
model=uipanel(figure1,'Units','normalized',...
        'Position',[0 0.7 1 0.3]);%left bottom width height
set(model,'Title','Model');

%Create axes
%Create axes

ax = axes('Parent',figure1,'OuterPosition',[0 0 1 0.7]);%left bottom width height
set(ax,'ButtonDownFcn',@replot);
%
replot(ax,1);


    function PlotModel(Wp,Wm)
        delete(findobj(model,'type','axes'))
        ImTrans(Wp,Wm,model);
    end

    function PlotLaplace(Wp,Wm,sm)
        A=SNRlaplace(srange,Wp,Wm,fp,w);
        eh=loglog(srange,Aenv,'g','LineWidth',2,'Parent',ax);
        hold(ax,'on');
        yl=[min(Aenv) max(Aenv)];
        xlim(ax,[min(srange) max(srange)]);
        ylim(ax,yl);
        xlabel(ax,'Laplace parameter, s');
        ylabel(ax,'Laplace transform, A(s)');
        lh=line([sm;sm],yl','Color','r','LineStyle','--','LineWidth',1.5,'Parent',ax);
        ah=loglog(srange,A,'b','LineWidth',1.5,'Parent',ax);
        hold(ax,'off');
        set(eh,'ButtonDownFcn',@replot);
        set(lh,'ButtonDownFcn',@replot);
        set(ah,'ButtonDownFcn',@replot);
        set(ax,'ButtonDownFcn',@replot);
        legend([eh ah lh],{'Envelope','Model','Max here'},'Location','best');
    end

    function [Wp,Wm]=FindOptModel(sm)
        [Wp,Wm]=FindOptL(sm,length(w),10,'DispExit',false);
    end

    function replot(~,~)
        sm=ginput(1);
        sm(2)=[];
        [Wp,Wm]=FindOptModel(sm);
        PlotModel(Wp,Wm);
        PlotLaplace(Wp,Wm,sm);
    end

end

