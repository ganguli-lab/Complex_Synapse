function [ newmodel ] = RJrunGUI( truemodel )
%BWUPDATEGUI Summary of this function goes here
%   Detailed explanation goes here



S.num_ch=7;
S.num_t=100;
S.num_BW=50;
S.fp=0.5;
chunks=[S.num_t*(0:S.num_ch-1)+1;S.num_t*(1:S.num_ch)];
AxFontSize=12;
BtFontSize=16;
cmapname='Hot';

%-------------------------------------------------------------------------
%Create figure
figure1 = figure('Units','normalized','OuterPosition',[0 0.03 1 0.97]);

%Create panels
model_true=uipanel(figure1,'Units','normalized',...
        'Position',[0 0 0.2 1]);%left bottom width height
set(model_true,'Title','True model');
%
model_est=uipanel(figure1,'Units','normalized',...
        'Position',[0.8 0 0.2 1]);%left bottom width height
set(model_est,'Title','Estimated model');
%
controls=uipanel(figure1,'Units','normalized',...
        'Position',[0.2 0.9 0.6 0.1]);%left bottom width height
set(controls,'Title','Controls');
%
simulation=uipanel(figure1,'Units','normalized',...
        'Position',[0.2 0.45 0.6 0.45]);%left bottom width height
set(simulation,'Title','Simulation');
%
metrics=uipanel(figure1,'Units','normalized',...
        'Position',[0.2 0 0.6 0.45]);%left bottom width height
set(metrics,'Title','Metrics');

%Create axes
ax_true(1) = axes('Parent',model_true,'OuterPosition',[0 0.6 1 0.4]);%left bottom width height
ax_true(2) = axes('Parent',model_true,'OuterPosition',[0 0.2 1 0.4]);%left bottom width height
ax_true(3) = axes('Parent',model_true,'OuterPosition',[0 0 1 0.2]);%left bottom width height

%
ax_est(1) = axes('Parent',model_est,'OuterPosition',[0 0.6 1 0.4]);%left bottom width height
ax_est(2) = axes('Parent',model_est,'OuterPosition',[0 0.2 1 0.4]);%left bottom width height
ax_est(3) = axes('Parent',model_est,'OuterPosition',[0 0 1 0.2]);%left bottom width height
%
potdepwt = axes('Parent',simulation,'OuterPosition',[0 0.8 0.9 0.2]);%left bottom width height
statepr  = axes('Parent',simulation,'OuterPosition',[0 0 0.95 0.8]);%left bottom width height
%
mets_true(1)=axes('Parent',metrics,'OuterPosition',[0 0.5 1 0.5]);%left bottom width height
mets_true(2)=axes('Parent',metrics,'OuterPosition',[0 0.5 1 0.5],'Color','none','YAxisLocation','right');%left bottom width height
mets_est(1)=axes('Parent',metrics,'OuterPosition',[0 0 1 0.5]);%left bottom width height
mets_est(2)=axes('Parent',metrics,'OuterPosition',[0 0 1 0.5],'Color','none','YAxisLocation','right');%left bottom width height

%where buttons and edit-boxes will be
edpos={0.0,0.0,0.5,1};%left bottom width height
btpos={0.5,0.0,0.5,1};%left bottom width height

pdleg=axes('Parent',simulation,'Position',[0.9 0.93 0.1 0.04]);
imagesc([1 -1],'Parent',pdleg);
set(pdleg,'XAxisLocation','top','XTick',[1 2],'YTick',[],'XTickLabel',{'pot','dep'});
wtleg=axes('Parent',simulation,'Position',[0.9 0.89 0.1 0.04]);
imagesc([1 -1],'Parent',wtleg);
set(wtleg,'XAxisLocation','bottom','XTick',[1 2],'YTick',[],'XTickLabel',{'strong','weak'});

colormap(pdleg,'jet')
freezeColors(pdleg);
colormap(wtleg,'jet')
freezeColors(wtleg);


colormap(potdepwt,'jet');
freezeColors(potdepwt);
colormap(statepr,cmapname);

%-------------------------------------------------------------------------
truemodel=truemodel.Sort(S.fp);
stpr=[];
Id=eye(length(truemodel.w));
simobj=SynapsePlastSeq;
newmodel=SynapseIdModel;
InitRand;
PlotModel(truemodel,ax_true);
prevmodel=newmodel;
metvals=zeros(S.num_BW,2,4);
%-------------------------------------------------------------------------
%Make edit boxes
Sprops=fieldnames(S);
for i =1:length(Sprops)
    MakeEditBox(i,Sprops{i});
end
numbt=3;
MakeButton(1,'Simulate',@Simulate);
MakeButton(2,'Initialise',@InitRand);
MakeButton(3,'Update',@Update);

% nummet=4;
% MakeTextBox(2,1,'True');
% MakeTextBox(3,1,'Est');
% MakeTextBox(1,2,'Likelihood');
% MakeTextBox(1,3,'KL: pot');
% MakeTextBox(1,4,'KL: dep');
% MakeTextBox(1,5,'KL: init');
% txh=zeros(2,nummet);
% for i=1:2
%     for j=1:nummet
%         txh(i,j)=MakeTextBox(i+1,j+1,'');
%     end
% end
%-------------------------------------------------------------------------
% UpdateMets;
%-------------------------------------------------------------------------

%  Callbacks for Gui

    %callback for edit-boxes
    function edb_callback(source,~,varname)
        S.(varname)=(str2double(get(source,'String')));
%         S.(varname)=str2double(get(source,'String'));
%         MakePlot;
    end
%-------------------------------------------------------------------------

%  Utility functions

 
%     function pos=CalcPosVert(which,maxwhich,left,bottom,width,height)
%     %Calculate position for panel in vertical grid of size [left bottom width height]
%         pos=[left bottom+height/maxwhich*(maxwhich-which) width height/maxwhich];
%     end

    function pos=CalcPosHorz(which,maxwhich,left,bottom,width,height)
     %Calculate position for panel in horizontal grid of size [left bottom width height]
       pos=[left+width/maxwhich*(which-1) bottom width/maxwhich height];
    end

    function MakeEditBox(which,varname)
        %control
        uicontrol(controls,'Style','edit',...
            'String',num2str(S.(varname)),...
            'Max',1,'Min',0,...
            'Units','normalized','FontSize',BtFontSize,...
            'Position',CalcPosHorz(2*which,2*length(Sprops),edpos{:}),...
            'Callback',{@edb_callback,varname});
        %labels
        uicontrol(controls,'Style','text',...
            'String',varname,...
            'Units','normalized','FontSize',BtFontSize,...
            'Position',CalcPosHorz(2*which-1,2*length(Sprops),edpos{:}));
    end%function MakeEditBox

    function MakeButton(which,varname,func)
        %control
        uicontrol(controls,'Style','pushbutton',...
                        'String',varname,...
                        'Units','normalized','FontSize',BtFontSize,...
                        'Position',CalcPosHorz(which,numbt,btpos{:}),...
                        'Callback',func);
    end%function MakeEditBox

    function PlotSim
        wt=[3-2*simobj.potdep 1; 2*simobj.readouts-3];
        cla(potdepwt);
        imagesc(wt,'Parent',potdepwt);
        set(potdepwt,'YTick',[1 2],'YTickLabel',{'pot/dep','weight'});
        xlabel(potdepwt,'Time','FontSize',AxFontSize);
        colormap(potdepwt,'jet');
        freezeColors(potdepwt);
        colormap(statepr,cmapname);
%         cla(statepr);
%         plot(st,'g','LineWidth',3,'Parent',statepr);
%         xlabel(statepr,'Time','FontSize',AxFontSize);
%         ylabel(statepr,'State','FontSize',AxFontSize);
    end

    function PlotStatePr
        cla(statepr);
        imagesc(stpr,'Parent',statepr);
        set(statepr,'CLim',[0 1]);
        cb=colorbar('peer',statepr);
        cblabel(cb,'Probability','FontSize',AxFontSize,'Rotation',270,'VerticalAlignment','bottom');
        hold(statepr,'on');
        plot(simobj.stateseq,'g','LineWidth',3,'Parent',statepr);
        xlabel(statepr,'Time','FontSize',AxFontSize);
        ylabel(statepr,'State','FontSize',AxFontSize);
        hold(statepr,'off');
    end

    function PlotModel(modelobj,ax)
        cla(ax(1));
        cla(ax(2));
        cla(ax(3));
        modelobj.image(ax(3),ax(1:2),'FontSize',AxFontSize);
        cb=colorbar('peer',ax(3),'location','SouthOutside');
        cblabel(cb,'Probability','FontSize',AxFontSize);
    end

    function Simulate(~,~)
        chunks=[S.num_t*(0:S.num_ch-1)+1;S.num_t*(1:S.num_ch)];
        simobj=truemodel.Simulate(S.fp,rand(2,S.num_t*S.num_t));
        PlotSim;
        PlotModel(truemodel,ax_true);
        %InitRand;
        %UpdateMets;
        stpr=StateProbs(newmodel,simobj);
        PlotStatePr;
        %InitRand;
    end

    function InitRand(~,~)
        n=length(truemodel.w);
        M_new={RandTrans(n)+Id,RandTrans(n)+Id};
        init_new=RandTrans(n)+Id;
        init_new=init_new(1,:);
        newmodel=SynapseIdModel(truemodel,'M',M_new,'Initial',init_new);
        newmodel=newmodel.Sort(S.fp);
        PlotModel(newmodel,ax_est);
        if ~isempty(simobj.potdep)
            stpr=StateProbs(newmodel,simobj);
            PlotStatePr;
        end
    end

    function mets=CalcMets
        if ~isempty(simobj.potdep)
            mets(1,1)=ChunkedHMMloglike(chunks,truemodel,simobj);
            mets(2,1)=ChunkedHMMloglike(chunks,newmodel,simobj);
        else
            mets(1:2,1)=NaN(2,1);
        end
        [mets(1,4),mets(1,2:3)]=truemodel.KLdivs(newmodel);
        [mets(2,4),mets(2,2:3)]=prevmodel.KLdivs(newmodel);
    end

    function PlotMets(upno)
        cla(mets_true(1));
        cla(mets_true(2));
        hl=plot(mets_true(1),(1:upno)',squeeze(metvals(1:upno,:,1)),'k');
        title(mets_true(1),'Comparison with true model','FontSize',AxFontSize);
        xlabel(mets_true(1),'Update #','FontSize',AxFontSize);
        ylabel(mets_true(1),'Log Likelihood','FontSize',AxFontSize);
        hkl=plot(mets_true(2),(1:upno)',squeeze(metvals(1:upno,1,2:4)));
        ylabel(mets_true(2),'KL divergence','FontSize',AxFontSize);
        FixDoublePlots(mets_true);
        set(hl(1),'LineStyle','--','LineWidth',2);
        legend(hl,{'loglikelihood(true)','loglikelihood(est)'},'location','west');
        legend(hkl,{'M^{pot}','M^{dep}','initial'},'location','east');
        cla(mets_est(1));
        cla(mets_est(2));
        hl=plot(mets_est(1),(1:upno-1)',squeeze(diff(metvals(1:upno,2,1))),'k');
        title(mets_est(1),'Comparison with previous estimate','FontSize',AxFontSize);
        xlabel(mets_est(1),'Update #','FontSize',AxFontSize);
        ylabel(mets_est(1),'\Delta Log Likelihood','FontSize',AxFontSize);
        hkl=plot(mets_est(2),(1:upno)',squeeze(metvals(1:upno,2,2:4)));
        ylabel(mets_est(2),'KL divergence','FontSize',AxFontSize);
        FixDoublePlots(mets_est);
        if upno>1
            legend(hl,{'\Delta loglikelihood'},'location','west');
        end
        legend(hkl,{'M^{pot}','M^{dep}','initial'},'location','east');
    end

    function FixDoublePlots(ax)
        xlim(ax(1),[0 S.num_BW]);
        xlim(ax(2),get(ax(1),'XLim'));
%         ylimits=get(ax(2),'YLim');
%         ytix=length(get(ax(1),'YTick'))-1;
%         yinc=(ylimits(2)-ylimits(1))/ytix;
        set(ax(1),'box','off');
        set(ax(2),'Color','none','YAxisLocation','right','XTick',[],'box','off','Position',get(ax(1),'Position'));
%         set(ax(2),'Color','none','YAxisLocation','right','XTick',[],'YTick',ylimits(1):yinc:ylimits(2),'Position',get(ax(1),'Position'));
    end

    function UpdateMets(upno)
        metvals(upno,:,:)=CalcMets;
    end

    function UpdateStep(upno)
        prevmodel=newmodel;
        [newmodel,stpr]=RJupdate(newmodel,simobj,chunks);
        newmodel=newmodel.Sort(S.fp);
        UpdateMets(upno);
        PlotMets(upno);
        PlotStatePr;
        PlotModel(newmodel,ax_est);
        drawnow;
%         pause;
    end

    function Update(~,~)
       if ~isempty(simobj.potdep)
        for upno=1:S.num_BW
            UpdateStep(upno);
        end
        PlotMets(S.num_BW);
        PlotStatePr;
        PlotModel(newmodel,ax_est);
        drawnow;
       end
    end

end

