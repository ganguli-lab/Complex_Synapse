function [ M_new,init_new ] = BWupdateGUI( M,init,w )
%BWUPDATEGUI Summary of this function goes here
%   Detailed explanation goes here



S.num_t=100;
S.fp=0.5;
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
        'Position',[0.2 0.3 0.6 0.6]);%left bottom width height
set(simulation,'Title','Simulation');
%
metrics=uipanel(figure1,'Units','normalized',...
        'Position',[0.2 0 0.6 0.3]);%left bottom width height
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
potdepwt = axes('Parent',simulation,'OuterPosition',[0 0.8 1 0.2]);%left bottom width height
statepr  = axes('Parent',simulation,'OuterPosition',[0 0 1 0.8]);%left bottom width height

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
outProj={WeightProj(1,w),WeightProj(-1,w)};
[M,w,ix]=SortByWtEta(M,w,S.fp);
init=init(ix);
Id=eye(length(w));
st=[];
pd=[];
rd=[];
stpr=[];
M_new={};
init_new=[];
InitRand;
PlotModel(M,init,ax_true);
M_prev=M_new;
init_prev=init_new;
init_obs=init;

%-------------------------------------------------------------------------
%Make edit boxes
Sprops=fieldnames(S);
for i =1:length(Sprops)
    MakeEditBox(i,Sprops{i});
end
numbt=3;
MakeButton(1,'Simulate',@Simulate)
MakeButton(2,'Initialise',@InitRand)
MakeButton(3,'Update',@Update)

nummet=4;
MakeTextBox(2,1,'True');
MakeTextBox(3,1,'Est');
MakeTextBox(1,2,'Likelihood');
MakeTextBox(1,3,'KL: pot');
MakeTextBox(1,4,'KL: dep');
MakeTextBox(1,5,'KL: init');
txh=zeros(2,nummet);
for i=1:2
    for j=1:nummet
        txh(i,j)=MakeTextBox(i+1,j+1,'');
    end
end
%-------------------------------------------------------------------------
UpdateMets;
%-------------------------------------------------------------------------

%  Callbacks for PertGui

    %callback for edit-boxes
    function edb_callback(source,~,varname)
        S.(varname)=(str2double(get(source,'String')));
%         S.(varname)=str2double(get(source,'String'));
%         MakePlot;
    end
%-------------------------------------------------------------------------

%  Utility functions

 
    function pos=CalcPosVert(which,maxwhich,left,bottom,width,height)
    %Calculate position for panel in vertical grid of size [left bottom width height]
        pos=[left bottom+height/maxwhich*(maxwhich-which) width height/maxwhich];
    end

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

    function h=MakeTextBox(row,which,initval)
        pos=CalcPosVert(row,3,0,0,1,1);
        pos=num2cell(pos);
        %labels
        h=uicontrol(metrics,'Style','text',...
            'String',initval,...
            'Units','normalized','FontSize',BtFontSize,...
            'Position',CalcPosHorz(which,nummet+1,pos{:}));
    end%function MakeEditBox

    function PlotSim
        wt=[2*pd-1 1; w(st)'];
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
        plot(st,'g','LineWidth',3,'Parent',statepr);
        xlabel(statepr,'Time','FontSize',AxFontSize);
        ylabel(statepr,'State','FontSize',AxFontSize);
        hold(statepr,'off');
    end

    function PlotModel(M,init,ax)
        cla(ax(1));
        cla(ax(2));
        cla(ax(3));
        imagesc(M{1},'Parent',ax(1));
        title(ax(1),'M^{pot}','FontSize',AxFontSize);
        xlabel(ax(1),'To','FontSize',AxFontSize);
        ylabel(ax(1),'From','FontSize',AxFontSize);
        imagesc(M{2},'Parent',ax(2));
        title(ax(2),'M^{dep}','FontSize',AxFontSize);
        xlabel(ax(2),'To','FontSize',AxFontSize);
        ylabel(ax(2),'From','FontSize',AxFontSize);
        imagesc(init,'Parent',ax(3));
        title(ax(3),'Initial','FontSize',AxFontSize);
        xlabel(ax(3),'State','FontSize',AxFontSize);
        set(ax,'CLim',[0 1]);
        cb=colorbar('peer',ax(3),'location','SouthOutside');
        cblabel(cb,'Probability','FontSize',AxFontSize);
    end

    function Simulate(~,~)
        [st,pd]=StateSequence(init,M,rand(2,S.num_t),S.fp);
        rd=(3-w(st)')/2;
        PlotSim;
        init_obs=init*outProj{rd(1)};
        init_obs=init_obs/sum(init_obs);
        PlotModel(M,init_obs,ax_true);
        %InitRand;
        UpdateMets;
        stpr=StateProbs(rd,init_new,outProj,M_new,pd);
        PlotStatePr;
    end

    function InitRand(~,~)
        n=length(w);
        M_new={RandTrans(n)+Id,RandTrans(n)+Id};
        init_new=RandTrans(n)+Id;
        init_new=init_new(1,:);
        [M_new,~,ix]=SortByWtEta(M_new,w,S.fp);
        init_new=init_new(ix);
        PlotModel(M_new,init_new,ax_est);
        if ~isempty(rd)
            stpr=StateProbs(rd,init_new,outProj,M_new,pd);
            PlotStatePr;
        end
    end

    function mets=CalcMets
        if ~isempty(rd)
            mets(1,1)=HMMlike(rd,init_obs,outProj,M,pd);
            mets(2,1)=HMMlike(rd,init_new,outProj,M_new,pd);
        else
            mets(1:2,1)=NaN(2,1);
        end
        mets(1,2)=KLdiv(M{1},M_new{1});
        mets(1,3)=KLdiv(M{2},M_new{2});
        mets(1,4)=KLdiv(init_obs,init_new);
        mets(2,2)=KLdiv(M_prev{1},M_new{1});
        mets(2,3)=KLdiv(M_prev{2},M_new{2});
        mets(2,4)=KLdiv(init_prev,init_new);
    end

    function UpdateMets
        mets=CalcMets;
        for ii=1:2
            for jj=1:nummet
                set(txh(ii,jj),'String',num2str(mets(ii,jj)));
            end
        end
    end

    function Update(~,~)
       if ~isempty(rd)
        M_prev=M_new;
        init_prev=init_new;
        [M_new,init_new,stpr]=BWupdate(rd,init_prev,outProj,M_prev,pd);
        [M_new,~,ix]=SortByWtEta(M_new,w,S.fp);
        init_new=init_new(ix);
        UpdateMets;
        PlotStatePr;
        PlotModel(M_new,init_new,ax_est);
       end
    end

end

