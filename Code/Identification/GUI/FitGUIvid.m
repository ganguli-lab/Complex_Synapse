function [ newmodel,loglike ] = FitGUIvid( imwriter,truemodel,options,num_ch,num_t )
%BWUPDATEGUI Summary of this function goes here
%   Detailed explanation goes here


% S=struct('MaxIter',100,'TolFun',NaN,'TolX',1e-4,'TolFunChange',1,...
%     'Penalty',1,'num_ch',30,'num_t',50);
% T=struct('Algorithm','BW','Weighter','RJ','Penaliser','No','ModelDiff','KL');
AxFontSize=12;
% BtFontSize=12;
cmapname='Hot';
options=SynapseOptimset(options,'PlotFcn',@PlotFcn,'GroundTruth',truemodel);
%-------------------------------------------------------------------------
%Create figure
figure1 = figure('Units','pixels','OuterPosition',[0 64 2560 1536]);

%Create panels
model_true=uipanel(figure1,'Units','normalized',...
        'Position',[0 0 0.2 1]);%left bottom width height
set(model_true,'Title','True model');
%
model_est=uipanel(figure1,'Units','normalized',...
        'Position',[0.8 0 0.2 1]);%left bottom width height
set(model_est,'Title','Estimated model');
%
% controls=uipanel(figure1,'Units','normalized',...
%         'Position',[0.2 0.9 0.6 0.1]);%left bottom width height
% set(controls,'Title','Controls');
%
simulation=uipanel(figure1,'Units','normalized',...
        'Position',[0.2 0.45 0.6 0.55]);%left bottom width height
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
mets_true(2)=axes('Parent',metrics,'OuterPosition',[0 0.5 1 0.5],'YAxisLocation','right');%left bottom width height
mets_true(1)=axes('Parent',metrics,'OuterPosition',[0 0.5 1 0.5],'Color','none');%left bottom width height
mets_est(2)=axes('Parent',metrics,'OuterPosition',[0 0 1 0.5],'YAxisLocation','right');%left bottom width height
mets_est(1)=axes('Parent',metrics,'OuterPosition',[0 0 1 0.5],'Color','none');%left bottom width height


% %where buttons and edit-boxes will be
% edpos={0.0,0.5,0.5,0.5};%left bottom width height
% edtpos={0.0,0.0,0.5,0.5};%left bottom width height
% btpos={0.5,0.0,0.5,1};%left bottom width height

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
truemodel=truemodel.Sort;
simobj=SynapsePlastSeqSim(1,num_ch);
newmodel=SynapseIdModel;
stopbtn=false;
InitRand;
PlotModel(truemodel,ax_true);
metvals=zeros(options.MaxIter,2,4);
Simulate;
drawnow;
PrintFrame(0);
Update;
%-------------------------------------------------------------------------
% %Make edit boxes
% Sprops=fieldnames(S);
% for i =1:length(Sprops)
%     MakeEditBox(i,Sprops{i});
% end
% Tprops=fieldnames(T);
% for i =1:length(Tprops)
%     MakeEditBoxT(i,Tprops{i});
% end
% numbt=4;
% MakeButton(1,'Simulate',@Simulate);
% MakeButton(2,'Initialise',@InitRand);
% MakeButton(3,'Update',@Update);
% MakeButton(4,'Stop',@Stop);

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
%     function edb_callback(source,~,varname)
%         S.(varname)=(str2double(get(source,'String')));
% %         MakePlot;
%     end
%     function edbt_callback(source,~,varname)
%         T.(varname)=get(source,'String');
% %         MakePlot;
%     end
%-------------------------------------------------------------------------

%  Utility functions

%-------------------------------------------------------------------------
%Functions that build interface elements
 
%     function pos=CalcPosVert(which,maxwhich,left,bottom,width,height)
%     %Calculate position for panel in vertical grid of size [left bottom width height]
%         pos=[left bottom+height/maxwhich*(maxwhich-which) width height/maxwhich];
%     end
% 
%     function pos=CalcPosHorz(which,maxwhich,left,bottom,width,height)
%      %Calculate position for panel in horizontal grid of size [left bottom width height]
%        pos=[left+width/maxwhich*(which-1) bottom width/maxwhich height];
%     end
% 
%     function MakeEditBox(which,varname)
%         %control
%         uicontrol(controls,'Style','edit',...
%             'String',num2str(S.(varname)),...
%             'Max',1,'Min',0,...
%             'Units','normalized','FontSize',BtFontSize,...
%             'Position',CalcPosHorz(2*which,2*length(Sprops),edpos{:}),...
%             'Callback',{@edb_callback,varname});
%         %labels
%         uicontrol(controls,'Style','text',...
%             'String',varname,...
%             'Units','normalized','FontSize',BtFontSize,...
%             'Position',CalcPosHorz(2*which-1,2*length(Sprops),edpos{:}));
%     end%function MakeEditBox
% 
%     function MakeEditBoxT(which,varname)
%         %control
%         uicontrol(controls,'Style','edit',...
%             'String',T.(varname),...
%             'Max',1,'Min',0,...
%             'Units','normalized','FontSize',BtFontSize,...
%             'Position',CalcPosHorz(2*which,2*length(Tprops),edtpos{:}),...
%             'Callback',{@edbt_callback,varname});
%         %labels
%         uicontrol(controls,'Style','text',...
%             'String',varname,...
%             'Units','normalized','FontSize',BtFontSize,...
%             'Position',CalcPosHorz(2*which-1,2*length(Tprops),edtpos{:}));
%     end%function MakeEditBox
% 
%     function MakeButton(which,varname,func)
%         %control
%         uicontrol(controls,'Style','pushbutton',...
%                         'String',varname,...
%                         'Units','normalized','FontSize',BtFontSize,...
%                         'Position',CalcPosHorz(which,numbt,btpos{:}),...
%                         'Callback',func);
%     end%function MakeEditBox
% 
%-------------------------------------------------------------------------
%Callback functions for buttons

    function Simulate(~,~)
        simobj=SynapsePlastSeqSim(1,num_ch);
        for jj=1:num_ch
            simobj(jj)=truemodel.Simulate(rand(2,num_t));
        end
        PlotSim;
        PlotModel(truemodel,ax_true);
        %InitRand;
        %UpdateMets;
        stpr=StateProbs(newmodel,simobj);
        if iscell(stpr)
            PlotStatePr(stpr);
        else
            PlotStatePr({stpr});
        end
        %InitRand;
    end

    function InitRand(~,~)
        newmodel=truemodel.Randomise;
        newmodel=newmodel.Sort;
        PlotModel(newmodel,ax_est);
        if ~isempty(simobj(1).potdep)
            stpr=StateProbs(newmodel,simobj);
            if iscell(stpr)
                PlotStatePr(stpr);
            else
                PlotStatePr({stpr});
            end
        end
    end

%     function Stop(~,~)
%         stopbtn=true;
%     end

    function Update(~,~)
       if ~isempty(simobj(1).potdep)
%            options=SynapseOptimset(catstruct(S,T),'PlotFcn',@PlotFcn,'GroundTruth',truemodel);
           [newmodel,loglike,exitflag,output]=FitSynapse(simobj,newmodel,options);
           warndlg(['Exit flag: ' int2str(exitflag) '. ' output.message]);
%            stopbtn=false;
       end
    end

%-------------------------------------------------------------------------
%Functions that update plots

    function PlotSim
%         wt=[3-2*[simobj.potdep]; 2*[simobj.readouts]-3];
        wt=[3-2*simobj(1).potdep; 2*simobj(1).readouts-3];
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

    function PlotStatePr(stpr)
        cla(statepr);
%         imagesc([stpr{:}],'Parent',statepr);
        imagesc(stpr{1},'Parent',statepr);
        set(statepr,'CLim',[0 1]);
        cb=colorbar('peer',statepr);
        cblabel(cb,'Probability','FontSize',AxFontSize,'Rotation',270,'VerticalAlignment','bottom');
        hold(statepr,'on');
        plot([simobj.stateseq],'g','LineWidth',3,'Parent',statepr);
        xlabel(statepr,'Time','FontSize',AxFontSize);
        ylabel(statepr,'State','FontSize',AxFontSize);
        hold(statepr,'off');
    end

    function PlotModel(modelobj,ax)
        cla(ax(1));
        cla(ax(2));
        cla(ax(3));
        modelobj.image(ax(3),ax(1:2),{'FontSize',AxFontSize});
        cb=colorbar('peer',ax(3),'location','SouthOutside');
        cblabel(cb,'Probability','FontSize',AxFontSize);
    end

    function mets=CalcMets(optimvals)
        mets(1,1)=optimvals.truth.fval;
        mets(2,1)=optimvals.fval;
        mets(1,4)=optimvals.truth.dInitial;
        mets(1,2:3)=optimvals.truth.dM;
        mets(2,4)=optimvals.prev.dInitial;
        mets(2,2:3)=optimvals.prev.dM;
    end

    function PlotMets(upno)
        cla(mets_true(1));
        cla(mets_true(2));
        %
        hl=plot(mets_true(1),(1:upno)',squeeze(metvals(1:upno,:,1)),'k');
        title(mets_true(1),'Comparison with true model','FontSize',AxFontSize);
        ylabel(mets_true(1),'Log Likelihood','FontSize',AxFontSize);
        set(hl(1),'LineStyle','--','LineWidth',2);
        legend(hl,{'loglikelihood(true)','loglikelihood(est)'},'location','west');
        %
        hkl=plot(mets_true(2),(1:upno)',squeeze(metvals(1:upno,1,2:4)));
        ylabel(mets_true(2),[options.ModelDiff ' distance'],'FontSize',AxFontSize);
        legend(hkl,{'M^{pot}','M^{dep}','initial'},'location','east');
        %
        xlabel(mets_true(1),'Update #','FontSize',AxFontSize);
        FixDoublePlots(mets_true);
        %
        cla(mets_est(1));
        cla(mets_est(2));
        %
        hl=plot(mets_est(1),(1:upno-1)',squeeze(diff(metvals(1:upno,2,1))),'k');
        title(mets_est(1),'Comparison with previous estimate','FontSize',AxFontSize);
        ylabel(mets_est(1),'\Delta Log Likelihood','FontSize',AxFontSize);
        if upno>1
            legend(hl,{'\Delta loglikelihood'},'location','west');
        end
        %
        hkl=plot(mets_est(2),(1:upno)',squeeze(metvals(1:upno,2,2:4)));
        ylabel(mets_est(2),[options.ModelDiff ' distance'],'FontSize',AxFontSize);
        legend(hkl,{'M^{pot}','M^{dep}','initial'},'location','east');
        %
        xlabel(mets_est(1),'Update #','FontSize',AxFontSize);
        FixDoublePlots(mets_est);
    end

    function FixDoublePlots(ax)
        xlim(ax(1),[0 options.MaxIter]);
        xlim(ax(2),get(ax(1),'XLim'));
%         ylimits=get(ax(2),'YLim');
%         ytix=length(get(ax(1),'YTick'))-1;
%         yinc=(ylimits(2)-ylimits(1))/ytix;
        set(ax(1),'box','off','Color','none');
        set(ax(2),'YAxisLocation','right','XTick',[],'box','off','Position',get(ax(1),'Position'));
%         set(ax(2),'Color','none','YAxisLocation','right','XTick',[],'YTick',ylimits(1):yinc:ylimits(2),'Position',get(ax(1),'Position'));
    end

    function UpdateMets(optimvals)
        metvals(optimvals.iteration,:,:)=CalcMets(optimvals);
    end

%     function UpdateStep(upno)
%         prevmodel=newmodel;
%         [newmodel,stpr]=RJweight(newmodel,simobj);
%         newmodel=newmodel.Sort(S.fp);
%         UpdateMets(upno);
%         PlotMets(upno);
%         PlotStatePr(stpr);
%         PlotModel(newmodel,ax_est);
%         drawnow;
% %         pause;
%     end
% 
    function stop=PlotFcn(model,optimValues,state)
        stop=stopbtn;
        switch state
            case 'init'
                set(figure1,'Name','Running');
            case 'iter'
                PlotModel(model,ax_est);
                PlotStatePr(optimValues.stateProb);
                UpdateMets(optimValues);
                PlotMets(optimValues.iteration);
                drawnow;
                PrintFrame(optimValues.iteration);
            case 'done'
                set(figure1,'Name','Done');
            otherwise
        end
    end

    function PrintFrame(imnumber)
    f=getframe(figure1);
    [im,map] = frame2im(f);    %Return associated image data 
    if isempty(map)            %Truecolor system
      rgb = im;
    else                       %Indexed system
      rgb = ind2rgb(im,map);   %Convert image data
    end
        imwriter.writeFrame(rgb,imnumber);
    end

end

