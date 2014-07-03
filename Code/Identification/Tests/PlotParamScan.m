function [ h ] = PlotParamScan( scan_results, paramvec, paramname )
%PLOTPARAMSCAN Summary of this function goes here
%   Detailed explanation goes here

doSize=scan_results(1).prob_st(1);
doDist=scan_results(1).KL(1);

h.fig=figure;
if doSize
    h.axSize=subplot(1,1+2*doDist,1);
    hold(h.axSize,'all');
end
if doDist
    h.axKL=subplot(1,doSize+2,doSize+1);
    hold(h.axKL,'all');
    h.axLn=subplot(1,doSize+2,doSize+2);
    hold(h.axLn,'all');
end

legtext=cell(length(scan_results),1);

for i=1:length(scan_results)
    legtext{i}=[paramname '=' num2str(paramvec(i))];
    if doSize
        errorbar(scan_results(i).num_events,scan_results(i).prob_st(1,:),scan_results(i).prob_st(2,:),'Parent',h.axSize);
    end
    if doDist
        errorbar(scan_results(i).num_events,scan_results(i).KL(1,:),scan_results(i).KL(2,:),'Parent',h.axKL);
        errorbar(scan_results(i).num_events,scan_results(i).Ln(1,:),scan_results(i).Ln(2,:),'Parent',h.axLn);
    end
end
    
if doSize
    ylim(h.axSize,[0 1]);
    xlabel(h.axSize,'# events');
    ylabel(h.axSize,'prob # states correct');
    legend(h.axSize,legtext);
    hold(h.axSize,'off');
end
if doDist
    xlabel(h.axKL,'# events');
    ylabel(h.axKL,'KL distance');
    legend(h.axKL,legtext);
    hold(h.axKL,'off');
    xlabel(h.axLn,'# events');
    ylabel(h.axLn,'L^n distance');
    legend(h.axLn,legtext);
    hold(h.axLn,'off');
end


end

