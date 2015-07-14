function LaplaceConstrEnvGUI( env, AenvSingle, NumSynapse )
%LAPLACECONSTRENVGUI(env, AenvSingle, NumSynapse) GUI for
%displaying constrained Laplacian envelope 
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Formatting options

AxFontSize=20;
BtFontSize=16;
EnvLinewidth=3;
ModelLineWidth=2;
TimeLineWidth=2;
EqLineWidth=2;
Interpreter='tex';
LumpThresh=2e-2;
DegThresh=1e-2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Data


AenvSingle=sqrt(NumSynapse)*AenvSingle./env.tau;

snr_xlim=[env.tau(1) env.tau(end)];
snr_ylim=[min(sqrt(NumSynapse)*env.SNRbenv(end),sqrt(NumSynapse)*env.sc*env.Ac) 5*AenvSingle(1)];

Play=false;
frNumber=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure

figure1 = figure('Units','pixels','OuterPosition',[0 64 2560 1536]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Panels

%snr curves
snr_ph=uipanel(figure1,'Units','normalized',...
    'Position',[0 0.15 0.7 0.85]);%left bottom width height
set(snr_ph,'Title','SNR curves','FontSize',AxFontSize);
%    'BackgroundColor','w',...

%current model
model_ph=uipanel(figure1,'Units','normalized',...
    'Position',[0.7 0.15 0.3 0.85]);%left bottom width height
set(model_ph,'Title','Current model','FontSize',AxFontSize);

%slider for frameNumber
slider_ph=uipanel(figure1,'Units','normalized',...
    'Position',[0 0 1 0.15],'FontSize',AxFontSize);
set(slider_ph,'Title','Frame');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Axes

%snr curves
snr_ax = axes('Parent',snr_ph,'OuterPosition',[0 0 1 1]);%left bottom width height

%synaptic weight heatmap
model_ax(1) = axes('Parent',model_ph,...
    'OuterPosition',[0 0.85 0.85 0.15]);%left bottom width height
%potentiation heatmap
model_ax(2) = axes('Parent',model_ph,...
    'OuterPosition',[0 0.5 0.85 0.35]);%left bottom width height
%depression heatmap
model_ax(3) = axes('Parent',model_ph,...
    'OuterPosition',[0 0.15 0.85 0.35]);%left bottom width height
%equilibrium distribution heatmap
model_ax(4) = axes('Parent',model_ph,...
    'OuterPosition',[0 0 0.85 0.15]);%left bottom width height

%colorbar for heatmaps
wt_cbpos=[0.85 0.9 0.05 0.08];
pr_cbpos=[0.85 0.05 0.05 0.8];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create slider: frameNumber
%slider
slh4=uicontrol(slider_ph,'Style','slider',...
                'Max',length(env.tau),'Min',1,...
                'Value',frNumber,...
                'SliderStep',[1 10]/length(env.tau),...
                'Units','normalized',...
                'Position',[0.1 0.25 0.8 0.7],...
                'Callback',{@sl_callback4});
%slider labels
uicontrol(slider_ph,'Style','text',...
                'String',env.tau(1),...
                'Units','normalized',...
                'FontSize',BtFontSize,...
                'Position',[0.1 0.05 0.08 0.2]);
uicontrol(slider_ph,'Style','text',...
                'String',env.tau(end),...
                'Units','normalized',...
                'FontSize',BtFontSize,...
                'Position',[0.85 0.05 0.05 0.2]);
%Play button
pbh4=uicontrol(slider_ph,'Style','pushbutton',...
                'String','Play',...
                'Units','normalized',...
                'Position',[0.0 0.25 0.1 0.7],...
                'FontSize',BtFontSize,...
                'Callback',{@pb_callback4});
%display value of slider
ed4=uicontrol(slider_ph,'Style','edit',...
                'String',num2str(env.tau(frNumber)),...
                'Units','normalized',...
                'Position',[0.9 0.6 0.1 0.35],...
                'FontSize',BtFontSize,...
                'Callback',{@ed_callback4});
%display value of slider
ed2=uicontrol(slider_ph,'Style','edit',...
                'String',num2str(frNumber),...
                'Units','normalized',...
                'Position',[0.9 0.25 0.1 0.35],...
                'FontSize',BtFontSize,...
                'Callback',{@ed_callback2});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','MATLAB:handle_graphics:exceptions:SceneNode');

changeFrameNumber(frNumber);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting functions

    function PlotSNR(frameNumber)
        cla(snr_ax);
        
        loglog(env.tau,AenvSingle,'g',env.tau,sqrt(NumSynapse)*env.SNRbenv,'r','Parent',snr_ax,'LineWidth',EnvLinewidth);
        hold(snr_ax,'on');
        loglog(1/env.sc,sqrt(NumSynapse)*env.sc*env.Ac,'ro','Parent',snr_ax,'MarkerFaceColor','r','MarkerSize',10);
        
        if ~isempty(env.mats(frameNumber).snrb)
            loglog(env.tau,sqrt(NumSynapse)*env.mats(frameNumber).snrb,'b','Parent',snr_ax,'LineWidth',ModelLineWidth);
        end
        
        line(env.tau(frameNumber)*[1;1],snr_ylim','Color','k','Parent',snr_ax,'LineWidth',TimeLineWidth)
        
        xlim(snr_ax,snr_xlim);
        ylim(snr_ax,snr_ylim);
        
        set(snr_ax,'TickLabelInterpreter',Interpreter,'FontSize',AxFontSize)
        
        switch Interpreter
            case 'tex'
                xlabel(snr_ax,'Mean recall time, \tau','FontSize',AxFontSize);
                ylabel(snr_ax,'\textsf{Recognition performance}, $\overline{\mathsf{SNR}}\mathsf{(\tau)}$','Interpreter','latex','FontSize',AxFontSize);
            case 'latex'
                xlabel(snr_ax,'Mean recall time, $\mathsf{\tau}$','Interpreter','latex','FontSize',AxFontSize);
                ylabel(snr_ax,'Recognition performance, $\overline{\mathrm{SNR}}(\tau)$','Interpreter','latex','FontSize',AxFontSize);
        end
        title(snr_ax,'Numerical envelopes','Interpreter',Interpreter,'FontSize',AxFontSize);
        legend(snr_ax,{'Unconstrained envelope',...
            'Constrained envelope',...
            'Constraint',...
            'Current model'},...
            'Location','northeast','Interpreter',Interpreter,'FontSize',AxFontSize);
    end

   function PlotModel(frameNumber)
       
       if env.mats(frameNumber).modelobj.isvalid
           modelobj=env.mats(frameNumber).modelobj;
           modelobj.LumpThresh=LumpThresh;
           modelobj.DegThresh=DegThresh;
           pt=modelobj.FindLumps;
    %        if modelobj.TestLump(pt)
               modelobj=modelobj.Lumpify(pt);
    %        end
           M=modelobj.NumStates;
           Id=eye(M);

           imagesc(modelobj.w','Parent',model_ax(1),[-1 1]);
           imagesc(modelobj.Wp+Id,'Parent',model_ax(2),[0 1]);
           imagesc(modelobj.Wm+Id,'Parent',model_ax(3),[0 1]);
           imagesc(modelobj.EqProb,'Parent',model_ax(4),[0 1]);
              
           line([1;1]*(1.5:1:M-0.5),[0.5;1.5]*ones(1,M-1),'Parent',model_ax(1),'Color',[0.5 0.5 0.5],'LineWidth',EqLineWidth);
           line([1;1]*(1.5:1:M-0.5),[0.5;M+1.5]*ones(1,M-1),'Parent',model_ax(2),'Color',[0.5 0.5 0.5],'LineWidth',EqLineWidth);
           line([0.5;M+1.5]*ones(1,M-1),[1;1]*(1.5:1:M-0.5),'Parent',model_ax(2),'Color',[0.5 0.5 0.5],'LineWidth',EqLineWidth);
           line([1;1]*(1.5:1:M-0.5),[0.5;M+1.5]*ones(1,M-1),'Parent',model_ax(3),'Color',[0.5 0.5 0.5],'LineWidth',EqLineWidth);
           line([0.5;M+1.5]*ones(1,M-1),[1;1]*(1.5:1:M-0.5),'Parent',model_ax(3),'Color',[0.5 0.5 0.5],'LineWidth',EqLineWidth);
           line([1;1]*(1.5:1:M-0.5),[0.5;1.5]*ones(1,M-1),'Parent',model_ax(4),'Color',[0.5 0.5 0.5],'LineWidth',EqLineWidth);
       else
           imagesc(0,'Parent',model_ax(1),[-1 1]);
           imagesc(0,'Parent',model_ax(2),[0 1]);
           imagesc(0,'Parent',model_ax(3),[0 1]);
           imagesc(0,'Parent',model_ax(4),[0 1]);
       end
       
       title(model_ax(1),'Synaptic weight','Interpreter',Interpreter,'FontSize',AxFontSize);
       title(model_ax(2),'Potentiation transition probability','Interpreter',Interpreter,'FontSize',AxFontSize);
       title(model_ax(3),'Depression transition probability','Interpreter',Interpreter,'FontSize',AxFontSize);
       title(model_ax(4),'Equilibrium probability','Interpreter',Interpreter,'FontSize',AxFontSize);
       
       xlabel(model_ax(1),'State','Interpreter',Interpreter,'FontSize',AxFontSize);
       xlabel(model_ax(2),'To state','Interpreter',Interpreter,'FontSize',AxFontSize);
       xlabel(model_ax(3),'To state','Interpreter',Interpreter,'FontSize',AxFontSize);
       xlabel(model_ax(4),'State','Interpreter',Interpreter,'FontSize',AxFontSize);
       
       ylabel(model_ax(2),'From state','Interpreter',Interpreter,'FontSize',AxFontSize);
       ylabel(model_ax(3),'From state','Interpreter',Interpreter,'FontSize',AxFontSize);

       colormap(model_ax(1),'redblue');
       colormap(model_ax(2),'hot');
       colormap(model_ax(3),'hot');
       colormap(model_ax(4),'hot');
 
       ch=colorbar('peer',model_ax(1),'Location','manual','Position',wt_cbpos,'TickLabelInterpreter',Interpreter,'FontSize',AxFontSize);
       colorbarlabel(ch,'Weight','Interpreter',Interpreter,'FontSize',AxFontSize);       

       ch=colorbar('peer',model_ax(4),'Location','manual','Position',pr_cbpos,'TickLabelInterpreter',Interpreter,'FontSize',AxFontSize);
       colorbarlabel(ch,'Probability','Interpreter',Interpreter,'FontSize',AxFontSize);       

       set(model_ax(1),'YTickMode','manual','YTick',[],'TickLabelInterpreter',Interpreter,'FontSize',AxFontSize);
       set(model_ax(2),'TickLabelInterpreter',Interpreter,'FontSize',AxFontSize);
       set(model_ax(3),'TickLabelInterpreter',Interpreter,'FontSize',AxFontSize);
       set(model_ax(4),'YTickMode','manual','YTick',[],'TickLabelInterpreter',Interpreter,'FontSize',AxFontSize);
   end

    function changeFrameNumber(frameNumber)
        frNumber=frameNumber;
        set(ed4,'String',num2str(env.tau(frameNumber)));
        set(ed2,'String',num2str(frNumber));
        set(slh4,'Value',frameNumber);
        PlotSNR(frameNumber);
        PlotModel(frameNumber);
        drawnow;
    end

    function DoPlay
       while Play && frNumber < length(env.tau)
           changeFrameNumber(frNumber+1);
       end %while Play
       Play=true;
       pb_callback4(pbh4,0);
    end %function DoPlay

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callbacks



    function sl_callback4(source,~)
        changeFrameNumber(round(get(source,'Value')));
    end %sl1_callback

    function ed_callback4(source,~)
        time=str2double(get(source,'String'));
        frameNumber=find(env.tau>=time,1,'first');
        changeFrameNumber(frameNumber);
    end %sl1_callback

    function ed_callback2(source,~)
        frameNumber=str2double(get(source,'String'));
        frameNumber=max(1,frameNumber);
        frameNumber=min(length(env.tau),frameNumber);
        changeFrameNumber(frameNumber);
    end %sl1_callback

    function pb_callback4(source,~)
        Play=~Play;
        if Play
            set(source,'String','Pause');
            DoPlay;
        else
            set(source,'String','Play');
        end %if
    end %sl1_callback




end

