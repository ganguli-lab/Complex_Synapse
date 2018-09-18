function PrintFigs( obj,prefix )
%VORexperiment.PRINTFIGS(prefix) printing VOR comparisons between WT/KO
%with/without pretraining
%   prefix = string to prepend to eps file names

obj.LabFontSize=2*obj.LabFontSize;
obj.txFontSize=2*obj.txFontSize;
obj.EqFontSize=2*obj.EqFontSize;
obj.ProbFontSize=4*obj.ProbFontSize;


    fig=figure('PaperPositionMode','auto','Position',[60 60 1000 1000]);
    figs=figure('PaperPositionMode','auto','Position',[60 60 1000 1000]);
    figEq=figure('PaperPositionMode','auto','Position',[60 60 1000 500]);
    figEv(1)=figure('PaperPositionMode','auto','Position',[60 60 600 500]);
    figEv(2)=figure('PaperPositionMode','auto','Position',[60 60 600 500]);
    
    Parent=axes('Parent',fig);
%     axes(Parent);
    obj.PlotLearn('LineWidth',2,'Parent',Parent);
    Parent=axes('Parent',figs);
%     axes(Parent);
    obj.PlotLearnS('LineWidth',2,'Parent',Parent);

    obj.EqProbPlots(figEq,'LineWidth',2);
    obj.ProbEvols(figEv);
    

    print(fig,[prefix '_learn.eps'],'-depsc');
    print(figs,[prefix '_learnS.eps'],'-depsc');
    print(figEq,[prefix '_eq_WT.eps'],'-depsc');
    print(figEv(1),[prefix '_pr_WT_nopre.eps'],'-depsc');
    print(figEv(2),[prefix '_pr_WT_pre.eps'],'-depsc');
    close('all');


end

