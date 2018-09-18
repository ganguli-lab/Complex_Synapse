function ViewFigs( obj )
%VORexperiment.VIEWFIGS viewing VOR comparisons between WT/KO
%with/without pretraining
%   Detailed explanation goes here


    fig=figure('WindowStyle','docked','PaperPositionMode','auto');
    figs=figure('WindowStyle','docked','PaperPositionMode','auto');
    figEq=figure('WindowStyle','docked','PaperPositionMode','auto');
    figEv=figure('WindowStyle','docked','PaperPositionMode','auto');

    Parent=axes('Parent',fig);
    obj.PlotLearn('LineWidth',2,'Parent',Parent);
    Parent=axes('Parent',figs);
    obj.PlotLearnS('LineWidth',2,'Parent',Parent);

    obj.EqProbPlots(figEq,'LineWidth',2);
    obj.ProbEvols(figEv);


end

