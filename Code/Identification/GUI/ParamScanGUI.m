function ParamScanGUI( scan_results,  paramName1,paramVals1,  paramName2,paramVals2 )
%ParamScanGUI(scan_results, paramName1,paramVals1, paramName2,paramVals2)
%GUI for displaying results of ScanParams
%   scan_results = output of ParamScan
%   paramName1/2 = name of parameter scanned for rows/cols of scan_results
%       (must be a valid struct field name)
%   paramVals1/2 = vector of values for parameter scanned for rows/cols of scan_results

S=struct('Metric',1,paramName1,1,paramName2,2);
mets={'prob_st','Ln','KL'};
metLabels={'prob # states correct','L^n distance','KL distance'};
LineStyles={'-','--','-.',':'};
Markers={'p','s','o','^','none'};


%-------------------------------------------------------------------------
%Create figure
figure1 = figure('Units','normalized','OuterPosition',[0 0.2 0.5 0.5]);

%Create panels
MakeRadioBox('Metric',['#states';'L^n    ';'KL     '],[0 0.8 1 0.2]);
MakeRadioBox(paramName2,paramVals2,[0 0 0.5 0.2]);
MakeRadioBox(paramName1,paramVals1,[0.5 0 0.5 0.2]);


plot_panel=uipanel(figure1,'Units','normalized',...
        'Position',[0 0.2 1 0.6]);%left bottom width height
%

%Create axes
ax(1) = axes('Parent',plot_panel,'FontSize',16,'OuterPosition',[0 0 0.5 1]);%left bottom width height
ax(2) = axes('Parent',plot_panel,'FontSize',16,'OuterPosition',[0.5 0 0.5 1]);%left bottom width height


MakePlots;



%-------------------------------------------------------------------------

%  Callbacks for Gui

    %callback for edit-boxes
    function sel_callback(~,eventdata,varname)
        S.(varname)=get(eventdata.NewValue,'UserData');
        MakePlots;
    end
%-------------------------------------------------------------------------

%  Utility functions

%-------------------------------------------------------------------------
%Functions that build interface elements
 
%     function pos=CalcPosVert(which,maxwhich,left,bottom,width,height)
%     %Calculate position for panel in vertical grid of size [left bottom width height]
%         pos=[left bottom+height/maxwhich*(maxwhich-which) width height/maxwhich];
%     end

    function pos=CalcPosHorz(which,maxwhich,left,bottom,width,height)
     %Calculate position for panel in horizontal grid of size [left bottom width height]
       pos=[left+width/maxwhich*(which-1) bottom width/maxwhich height];
    end

    function MakeRadioBox(varname,varvals,pos)
        panel=uipanel(figure1,'Title',varname,'Units','normalized','Position',pos,'FontSize',16);
        btgroup=uibuttongroup('Parent',panel,'Units','normalized','Visible','off',...
            'Position',[0 0 1 1],'SelectionChangeFcn',{@sel_callback,varname});
        %control
        if ischar(varvals)
            varstr=varvals;
        else
            varstr=num2str(varvals');
        end
        for i=1:size(varstr,1)
            uicontrol(btgroup,'Style','radiobutton','String',varstr(i,:),'UserData',i,...
                'Units','normalized','Position',CalcPosHorz(i,size(varstr,1),0,0,1,1),...
                'FontSize',20);
        end
        set(btgroup,'Visible','on');
    end%function MakeEditBox


%-------------------------------------------------------------------------
%Functions that update plots

    function MakePlots
        MakeOnePlot(ax(1),scan_results(:,S.(paramName2)),paramName1,paramVals1,paramName2,paramVals2(S.(paramName2)))
        MakeOnePlot(ax(2),scan_results(S.(paramName1),:),paramName2,paramVals2,paramName1,paramVals1(S.(paramName1)))
    end

    function MakeOnePlot(axh,data,varname,varvals,othervarname,othervarval)
        cla(axh,'reset');
        set(axh,'FontSize',16)
        hold(axh,'all');
        legtext=cell(1,length(varvals));
        for i=1:length(varvals)
            legtext{i}=[varname '=' num2str(varvals(i))];
            errorbar(data(i).num_events,data(i).(mets{S.Metric})(1,:),data(i).(mets{S.Metric})(2,:),...
                'Parent',axh,'LineWidth',2,...
                'LineStyle',LineStyles{mod(i-1,length(LineStyles))+1},...
                'Marker',Markers{mod(i-1,length(Markers))+1});
        end
    
        switch S.Metric
            case 1
                ylim(axh,[-0.05 1.05]);
            case 2
                ylim(axh,[0 0.9]);
            case 3
                ylim(axh,[0 0.2]);
        end
        title(axh,[othervarname '=' num2str(othervarval)]);
        xlabel(axh,'# events');
        ylabel(axh,metLabels{S.Metric});
        legend(axh,legtext,'Location','Best');
        hold(axh,'off');
    end%function MakeOnePlot




 end%function ParamScanGUI