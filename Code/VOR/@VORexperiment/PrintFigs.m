function PrintFigs( obj,prefix )
%PRINTFIGS Summary of this function goes here
%   Detailed explanation goes here

obj.LabFontSize=2*obj.LabFontSize;
obj.txFontSize=2*obj.txFontSize;
obj.EqFontSize=2*obj.EqFontSize;
obj.ProbFontSize=4*obj.ProbFontSize;


    fig=figure('PaperPositionMode','auto','Position',[60 60 1000 1000]);
    figs=figure('PaperPositionMode','auto','Position',[60 60 1000 1000]);
    figEq(1)=figure('PaperPositionMode','auto','Position',[60 60 1000 500]);
    figEq(2)=figure('PaperPositionMode','auto','Position',[60 60 1000 500]);
    figEv(1)=figure('PaperPositionMode','auto','Position',[60 60 600 500]);
    figEv(2)=figure('PaperPositionMode','auto','Position',[60 60 600 500]);
    figEv(3)=figure('PaperPositionMode','auto','Position',[60 60 600 500]);
    figEv(4)=figure('PaperPositionMode','auto','Position',[60 60 600 500]);
    
    Parent=axes('Parent',fig);
    obj.PlotLearn('LineWidth',2,'Parent',Parent);
    Parent=axes('Parent',figs);
    obj.PlotLearnS('LineWidth',2,'Parent',Parent);

    obj.EqProbPlots(figEq,'LineWidth',2);
    obj.ProbEvols(figEv);
    

    print(fig,[prefix '_learn.eps'],'-depsc');
    print(figs,[prefix '_learnS.eps'],'-depsc');
    print(figEq(1),[prefix '_eq_WT.eps'],'-depsc');
    print(figEq(2),[prefix '_eq_KO.eps'],'-depsc');
    print(figEv(1),[prefix '_pr_WT_nopre.eps'],'-depsc');
    print(figEv(2),[prefix '_pr_WT_pre.eps'],'-depsc');
    print(figEv(3),[prefix '_pr_KO_nopre.eps'],'-depsc');
    print(figEv(4),[prefix '_pr_KO_pre.eps'],'-depsc');
    close('all');


end

