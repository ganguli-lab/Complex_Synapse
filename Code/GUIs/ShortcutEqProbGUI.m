function ShortcutEqProbGUI(varargin)

% %initial params
% 
% varargin=assignApplicable(varargin);
% 
% 
% %check validity of params
% assert(isscalar(nstates));
% assert(isint(nstates));



%Initialise data (make them global)
S.nstates=4;
S.start=2;
S.len=2;
S.fp=1/2;
S.transprob=1;
S.pertval=0;
numsl=3;
numed=3;

%-------------------------------------------------------------------------
%Create figure
figure1 = figure('Units','normalized');

%Create axes 1
axes1 = axes('Parent',figure1,'OuterPosition',[0 0.2 1 0.8]);%left bottom width height


%-------------------------------------------------------------------------
%Make dummy handles
phs=zeros(1,numsl);
slh=phs;
sth=phs;
phe=zeros(1,numed);
%-------------------------------------------------------------------------
%Create param panels
MakeSlider(1,'pertval');
MakeSlider(2,'transprob');
MakeSlider(3,'fp');
MakeEditBox(1,'nstates')
MakeEditBox(2,'start')
MakeEditBox(3,'len')
%-------------------------------------------------------------------------
%Make initial plots
MakePlot;

%-------------------------------------------------------------------------

%  Callbacks for PertGui

    function sl_callback(source,~,which,varname)
        S.(varname)=get(source,'Value');
        set(sth(which),'String',num2str(S.(varname)));
        MakePlot;
    end 

    function pb_callback(~,~,which,varname)
        set(slh(which),'Value',0.5);
        sl_callback(slh(which),0,varname);
    end

    function edb_callback(source,~,varname)
        S.(varname)=floor(str2double(get(source,'String')));
        MakePlot;
    end


%-------------------------------------------------------------------------

%  Utility functions

    function MakeSlider(which,varname)
        phs(which)=uipanel(figure1,'Units','normalized',...
                'Position',[0.05 0.2/numsl*(numsl-which) 0.8 0.2/numsl]);%left bottom width height
        set(phs(which),'Title',varname);
        %slider
        slh(which)=uicontrol(phs(which),'Style','slider',...
                        'Max',1,'Min',0,'Value',S.(varname),...
                        'SliderStep',[0.02 0.2],...
                        'Units','normalized',...
                        'Position',[0.1 0.25 0.8 0.7],...
                        'Callback',{@sl_callback,which,varname});
        %slider labels
        uicontrol(phs(which),'Style','text',...
                        'String','0',...
                        'Units','normalized',...
                        'Position',[0.1 0.05 0.08 0.2]);
        uicontrol(phs(which),'Style','text',...
                        'String','0.5',...
                        'Units','normalized',...
                        'Position',[0.475 0.05 0.05 0.2]);
        uicontrol(phs(which),'Style','text',...
                        'String','1',...
                        'Units','normalized',...
                        'Position',[0.85 0.05 0.05 0.2]);
        %button to zero slider
        uicontrol(phs(which),'Style','pushbutton',...
                        'String','1/2',...
                        'Units','normalized',...
                        'Position',[0.0 0.25 0.1 0.7],...
                        'Callback',{@pb_callback,which,varname});
        %display value of slider
        sth(which)=uicontrol(phs(which),'Style','text',...
                        'String',num2str(get(slh(which),'Value')),...
                        'Units','normalized',...
                        'Position',[0.9 0.25 0.1 0.7]);
    end%function MakeSlider

    function MakeEditBox(which,varname)
        phe(which)=uipanel(figure1,'Units','normalized',...
                'Position',[0.87 0.2/numed*(numed-which) 0.1 0.2/numed]);%left bottom width height
        set(phe(which),'Title',varname);
        %control
        uicontrol(phe(which),'Style','edit',...
            'String',num2str(S.(varname)),...
            'Max',1,'Min',0,...
            'Units','normalized',...
            'Position',[0.05 0.05 0.9 0.9],...
            'Callback',{@edb_callback,varname});
        %labels
%         uicontrol(phe(which),'Style','text',...
%             'String',title,...
%             'Units','normalized',...
%             'Position',[x y 1/6 1/2]);
    end%function MakeEditBox

    function MakePlot
        q=ones(1,S.nstates-1)*S.transprob;
%         [wp,wm]=MakeSMS(q);
%         W=S.fp*wp+(1-S.fp)*wm;
%         [pertp,pertm]=ShortcutPert(W,S.start,S.len);
%         pert=S.fp*pertp+(1-S.fp)*pertm;
%         eqp=EqProb(W+S.pertval*pert);
        [Wp,Wm]=CascadeMSinterp(q,q,S.fp,S.pertval);
        W=S.fp*Wp+(1-S.fp)*Wm;
        eqp=EqProb(W);
        bar(axes1,eqp);
    end%function MakePlot



end %function IntGui