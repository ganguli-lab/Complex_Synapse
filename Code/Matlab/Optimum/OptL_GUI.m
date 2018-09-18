function OptL_GUI( fp,w,Aenv,srange )
%OPTL_GUI Summary of this function goes here
%   Detailed explanation goes here

partitions=num2cell(1:length(w));

%-------------------------------------------------------------------------
%Create figure
% figure1 = figure('Units','normalized','OuterPosition',[0 0.03 1 0.97]);
figure1 = figure('Units','normalized','WindowStyle','docked');

%Create panels
model=uipanel(figure1,'Units','normalized',...
        'Position',[0.5 0.66 0.5 0.33]);%left bottom width height
set(model,'Title','Model');
lumptest=uipanel(figure1,'Units','normalized',...
        'Position',[0.5 0.33 0.5 0.33]);%left bottom width height
set(lumptest,'Title','Lumapbility');
lumpmodel=uipanel(figure1,'Units','normalized',...
        'Position',[0.5 0 0.5 0.33]);%left bottom width height
set(lumpmodel,'Title','Lumped model');

%Create axes

ax = axes('Parent',figure1,'OuterPosition',[0 0 0.5 1]);%left bottom width height
set(ax,'ButtonDownFcn',@replot);
%
replot(ax,1);


    function PlotModel(Wp,Wm,panel,varargin)
        delete(findobj(panel,'type','axes'))
        ih=ImTrans(Wp,Wm,panel,[],[],varargin{:});
        set(ih,'ButtonDownFcn',@Lump);
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
        PlotModel(Wp,Wm,model);
        TestLump(Wp,Wm);
        PlotLump(Wp,Wm);
        PlotLaplace(Wp,Wm,sm);
    end

    function Lump(~,~)
        part_str=inputdlg({'Separate partitions with ;'},'Enter partitions',1,{int2str(1:length(w))});
        partitions=Str2partitions(part_str{1});
        ih=findobj(model,'type','image');
        Wp=get(ih(4),'CData')-eye(length(w));
        Wm=get(ih(3),'CData')-eye(length(w));
        TestLump(Wp,Wm);
        PlotLump(Wp,Wm);
    end

    function TestLump(Wp,Wm)
        Wptest=LumpTest(Wp+eye(length(Wp)),partitions);
        Wmtest=LumpTest(Wm+eye(length(Wm)),partitions);
        PlotModel(Wptest-eye(size(Wptest)),Wmtest-eye(size(Wmtest)),lumptest,'CLim',[]);
    end
        
    function PlotLump(Wp,Wm)
        Wptest=LumpedMat(Wp,partitions);
        Wmtest=LumpedMat(Wm,partitions);
        PlotModel(Wptest,Wmtest,lumpmodel);
    end
        

end

