function ParamScanGUI( scan_results,  paramName1,paramVals1,  paramName2,paramVals2 )

S=struct('Metric',1,paramName1,1,paramName2,2);
mets={'prob_st','KL','Ln'};
metLabels={'prob # states correct','KL distance','L^n distance'};

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
        cla(axh);
        hold(axh,'all');
        legtext=cell(1,length(varvals));
        for i=1:length(varvals)
            legtext{i}=[varname '=' num2str(varvals(i))];
            errorbar(data(i).num_events,data(i).(mets{S.Metric})(1,:),data(i).(mets{S.Metric})(2,:),'Parent',axh,'LineWidth',2);
        end
    
        if S.Metric==1
            ylim(axh,[0 1]);
        end
        title(axh,[othervarname '=' num2str(othervarval)]);
        xlabel(axh,'# events');
        ylabel(axh,metLabels{S.Metric});
        legend(axh,legtext);
        hold(axh,'off');
    end%function MakeOnePlot




 end%function ParamScanGUI