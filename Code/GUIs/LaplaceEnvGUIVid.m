function LaplaceEnvGUIVid( imwriter, chains, NumSynapse, inds )
%LAPLACEENVGUIVID(imwriter,chains,NumSynapse,inds) GUI for displaying Laplacian
%envelope 
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Formatting options

AxFontSize=30;
% BtFontSize=16;
EnvLinewidth=3;
ModelLineWidth=2;
TimeLineWidth=2;
EqLineWidth=2;
Interpreter='tex';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Data

if ~exist('inds','var') || isempty(inds)
    inds=1:length(chains);
end

chains=chains(inds);
srange=[chains.s];
Qmat=reshape([chains.qv],[],length(chains))';
AenvNum=[chains.A];

tau=1./srange;
M=size(Qmat,2)/2+1;

AenvNum=sqrt(NumSynapse)*AenvNum./tau;
AenvProven=sqrt(NumSynapse)*(M-1)./(tau+(M-1));

snr_xlim=[tau(1) tau(end)];
snr_ylim=[AenvNum(end) 5*sqrt(NumSynapse)];

% Play=false;
% frNumber=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure

figure1 = figure('Units','pixels','Position',[0 64 2560 1008],'PaperPositionMode','auto');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Panels

%snr curves
snr_ph=uipanel(figure1,'Units','normalized',...
    'BackgroundColor','w',...
    'Position',[0 0 0.5 1]);%left bottom width height
% set(snr_ph,'Title','SNR curves','FontSize',AxFontSize);

%current model
model_ph=uipanel(figure1,'Units','normalized',...
    'BackgroundColor','w',...
    'Position',[0.5 0 0.5 1]);%left bottom width height
% set(model_ph,'Title','Current model','FontSize',AxFontSize);

% %slider for frameNumber
% slider_ph=uipanel(figure1,'Units','normalized',...
%     'Position',[0 0 1 0.15],'FontSize',AxFontSize);
% set(slider_ph,'Title','Frame');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Axes

%snr curves
snr_ax = axes('Parent',snr_ph,'OuterPosition',[0 0 1 1]);%left bottom width height

%potentiation bar
model_ax(1) = axes('Parent',model_ph,...
    'XTickMode','manual',...
    'YDir','reverse',...
    'Position',[0.05 0.55 0.75 0.35]);%left bottom width height
%equilibrium distribution heatmap
model_ax(2) = axes('Parent',model_ph,...
    'XTickMode','manual',...
    'Position',[0.05 0.45 0.75 0.1]);%left bottom width height
%depression bar
model_ax(3) = axes('Parent',model_ph,...
    'XTickMode','manual',...
    'Position',[0.05 0.1 0.75 0.35]);%left bottom width height

%colorbar for equilibrium distribution heatmap
cbpos=[0.85 0.05 0.05 0.9];


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Create slider: frameNumber
% %slider
% slh4=uicontrol(slider_ph,'Style','slider',...
%                 'Max',length(tau),'Min',1,...
%                 'Value',frNumber,...
%                 'SliderStep',[1 10]/length(tau),...
%                 'Units','normalized',...
%                 'Position',[0.1 0.25 0.8 0.7],...
%                 'Callback',{@sl_callback4});
% %slider labels
% uicontrol(slider_ph,'Style','text',...
%                 'String',tau(1),...
%                 'Units','normalized',...
%                 'FontSize',BtFontSize,...
%                 'Position',[0.1 0.05 0.08 0.2]);
% uicontrol(slider_ph,'Style','text',...
%                 'String',tau(end),...
%                 'Units','normalized',...
%                 'FontSize',BtFontSize,...
%                 'Position',[0.85 0.05 0.05 0.2]);
% %Play button
% pbh4=...
% uicontrol(slider_ph,'Style','pushbutton',...
%                 'String','Play',...
%                 'Units','normalized',...
%                 'Position',[0.0 0.25 0.1 0.7],...
%                 'FontSize',BtFontSize,...
%                 'Callback',{@pb_callback4});
% %display value of slider
% ed4=uicontrol(slider_ph,'Style','edit',...
%                 'String',num2str(tau(frNumber)),...
%                 'Units','normalized',...
%                 'Position',[0.9 0.25 0.1 0.7],...
%                 'FontSize',BtFontSize,...
%                 'Callback',{@ed_callback4});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','MATLAB:handle_graphics:exceptions:SceneNode');

% changeFrameNumber(frNumber);
% 
%    while frNumber < length(tau)
%        changeFrameNumber(frNumber+1);
%    end %while Play
   
for i=1:length(tau)
    changeFrameNumber(i);
end

close(figure1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting functions


    function PlotSNR(frameNumber)
        cla(snr_ax);
        
        loglog(tau,AenvProven,'g',tau,AenvNum,'r','Parent',snr_ax,'LineWidth',EnvLinewidth);
        hold(snr_ax,'on');
        
        loglog(tau,sqrt(NumSynapse)*chains(frameNumber).snrb,'b','Parent',snr_ax,'LineWidth',ModelLineWidth);
        
        line(tau(frameNumber)*[1;1],snr_ylim','Color','k','Parent',snr_ax,'LineWidth',TimeLineWidth)
        
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
        title(snr_ax,'Proven and numerical envelopes','Interpreter',Interpreter,'FontSize',AxFontSize);
        legend(snr_ax,{'Proven envelope',...
            'Numeric envelope',...
            'Current model'},...
            'Location','northeast','Interpreter',Interpreter,'FontSize',AxFontSize);
    end

   function PlotModel(frameNumber)
       pr=chains(frameNumber).modelobj.EqProb;
       qp=Qmat(frameNumber,1:M-1);
       qm=Qmat(frameNumber,M:end);
       [qp,qm]=ZeroIrrel(qp,qm);
       
       bar(model_ax(1),qp,'r');
       imagesc(pr,'Parent',model_ax(2),[0 0.5]);
       line([1;1]*(1.5:1:M-0.5),[0.5;1.5]*ones(1,M-1),'Parent',model_ax(2),'Color',[0.5 0.5 0.5],'LineWidth',EqLineWidth);
       bar(model_ax(3),qm,'b');
       
       xlim(model_ax(1),[0 M]);
       xlim(model_ax(2),[0.5 M+0.5]);
       xlim(model_ax(3),[0 M]);
       ylim(model_ax(1),[0 1]);
       ylim(model_ax(2),[0.5 1.5]);
       ylim(model_ax(3),[0 1]);
       
       title(model_ax(1),'Potentiation transition probability','Interpreter',Interpreter,'FontSize',AxFontSize);
       th=title(model_ax(3),'Depression transition probability','Interpreter',Interpreter,'FontSize',AxFontSize);
       tpos=get(th,'Position');
       
       colormap(model_ax(2),'hot');
       ch=colorbar(model_ax(2),'Location','manual','TickLabelInterpreter',Interpreter,'Position',cbpos,'FontSize',AxFontSize);
       colorbarlabel(ch,'Equilibrium distribution','Interpreter',Interpreter,'FontSize',AxFontSize);
       
       set(model_ax(1),'XTickMode','manual','XTick',[],'TickLabelInterpreter',Interpreter,'FontSize',AxFontSize);
       set(model_ax(2),'XTickMode','manual','XTick',[],'YTickMode','manual','YTick',[],'TickLabelInterpreter',Interpreter,'FontSize',AxFontSize);
       set(model_ax(3),'XTickMode','manual','XTick',[],'YDir','reverse','TickLabelInterpreter',Interpreter,'FontSize',AxFontSize);
       set(th,'Position',tpos,'VerticalAlignment','top');
%        title(model_ax(3),'Depression transition probability','FontSize',AxFontSize);
       
   end

    function changeFrameNumber(frameNumber)
%         frNumber=frameNumber;
%         set(ed4,'String',num2str(tau(frameNumber)));
%         set(slh4,'Value',frameNumber);
        PlotSNR(frameNumber);
        PlotModel(frameNumber);
        drawnow;
        PrintFrame(frameNumber);
    end

%     function DoPlay
%        while Play && frNumber < length(tau)
%            changeFrameNumber(frNumber+1);
%        end %while Play
%        Play=true;
%        pb_callback4(pbh4);
%     end %function DoPlay
% 
    function [newqp,newqm]=ZeroIrrel(qp,qm)
        newqp=qp;
        newqm=qm;
        [qmin,ix]=min(qp);
        if qmin < 1.3e-3 && ix > M/2
            newqp(ix+1:end)=qmin;
            newqm(ix:end)=qmin;
        end
        [qmin,ix]=min(qm);
        if qmin < 1.3e-3 && ix < M/2
            newqp(1:ix)=qmin;
            newqm(1:ix-1)=qmin;
        end
    end

    function PrintFrame(frameNumber)
        imwriter.writeFig(figure1,frameNumber);
    end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callbacks



%     function sl_callback4(source,~)
%         changeFrameNumber(round(get(source,'Value')));
%     end %sl1_callback
% 
%     function ed_callback4(source,~)
%         time=str2double(get(source,'String'));
%         frameNumber=find(tau>=time,1,'first');
%         changeFrameNumber(frameNumber);
%     end %sl1_callback
% 
%     function pb_callback4(source,~)
%         Play=~Play;
%         if Play
%             set(source,'String','Pause');
%             DoPlay;
%         else
%             set(source,'String','Play');
%         end %if
%     end %sl1_callback
% 



end

