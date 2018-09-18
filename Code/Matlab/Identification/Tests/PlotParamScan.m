function [ figure1 ] = PlotParamScan( scan_results, paramName1,paramVals1,  paramName2,paramVals2, whichmet,whichval1,whichval2 )
%fh=PLOTPARAMSCAN(scan_results,paramName1,paramVals1,paramName2,paramVals2,whichmet,whichval1,whichval2)
%display results of ScanParams 
%   fh = figure handle
%   scan_results = output of ParamScan
%   paramName1/2 = name of parameter scanned for rows/cols of scan_results
%       (must be a valid struct field name)
%   paramVals1/2 = vector of values for parameter scanned for rows/cols of scan_results
%   whichmet = (1,2,3) to plot (probability # states correct,KL distance,L^n distance)
%   whichval1/2 = index of paramVals1/2 to use for plot of scan of paramName2/1

mets={'prob_st','KL','Ln'};
metLabels={'probability # states correct','KL distance','L^n distance'};
LineStyles={'-','--','-.',':'};
Markers={'p','s','o','^','none'};


%-------------------------------------------------------------------------
%Create figure
figure1 = figure('Units','normalized','OuterPosition',[0 0.2 0.5 0.5]);

%

%Create axes
ax(1) = axes('Parent',figure1,'FontSize',16,'OuterPosition',[0 0 0.5 1]);%left bottom width height
ax(2) = axes('Parent',figure1,'FontSize',16,'OuterPosition',[0.5 0 0.5 1]);%left bottom width height


MakeOnePlot(ax(1),scan_results(:,whichval2),paramName1,paramVals1,paramName2,paramVals2(whichval2))
MakeOnePlot(ax(2),scan_results(whichval1,:),paramName2,paramVals2,paramName1,paramVals1(whichval1))




%-------------------------------------------------------------------------
%Functions that make plots


    function MakeOnePlot(axh,data,varname,varvals,othervarname,othervarval)
        cla(axh);
        hold(axh,'all');
        legtext=cell(1,length(varvals));
        for i=1:length(varvals)
            legtext{i}=[varname '=' num2str(varvals(i))];
            errorbar(data(i).num_events,data(i).(mets{whichmet})(1,:),data(i).(mets{whichmet})(2,:),...
                'Parent',axh,'LineWidth',2,...
                'LineStyle',LineStyles{mod(i-1,length(LineStyles))+1},...
                'Marker',Markers{mod(i-1,length(Markers))+1});
        end
    
        if whichmet==1
            ylim(axh,[-0.05 1.05]);
        end
        title(axh,[othervarname '=' num2str(othervarval)]);
        xlabel(axh,'# events');
        ylabel(axh,metLabels{whichmet});
        legend(axh,legtext,'Location','Best');
        hold(axh,'off');
    end%function MakeOnePlot



end

