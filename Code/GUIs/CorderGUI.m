function CorderGUI(varargin)

% %initial params
nstates=2;
% 
% varargin=assignApplicable(varargin);
% 
% 
% %check validity of params
% assert(isscalar(nstates));
% assert(isint(nstates));



%Initialise data (make them global)
fp=1/2;%last element=f+

%-------------------------------------------------------------------------
%Create figure
figure1 = figure('Units','normalized');

%Create axes 1
axes1 = axes('Parent',figure1,'OuterPosition',[0 0.2 1 0.8]);%left bottom width height


%-------------------------------------------------------------------------
%Make dummy handles
ph=[0,0];
slh=0;
sth=0;
%-------------------------------------------------------------------------
%Create param panels
MakeSlider(1);
MakeEditBox(2)
%-------------------------------------------------------------------------
%Make initial plots
MakePlot;

%-------------------------------------------------------------------------

%  Callbacks for PertGui

    function sl_callback(source,~)
        fp=get(source,'Value');
        set(sth,'String',num2str(fp));
        MakePlot;
    end 

    function pb_callback(~,~)
        set(slh,'Value',0.5);
        sl_callback(slh,0);
    end

    function edb_callback(source,~)
        nstates=floor(str2double(get(source,'String')));
        MakePlot;
    end


%-------------------------------------------------------------------------

%  Utility functions

    function MakeSlider(which)
        ph(which)=uipanel(figure1,'Units','normalized',...
                'Position',[0.05 0.05 0.8 0.1]);%left bottom width height
        set(ph(which),'Title','Fraction of potentiation');
        %slider
        slh(which)=uicontrol(ph(which),'Style','slider',...
                        'Max',1,'Min',0,'Value',fp,...
                        'SliderStep',[0.02 0.2],...
                        'Units','normalized',...
                        'Position',[0.1 0.25 0.8 0.7],...
                        'Callback',{@sl_callback});
        %slider labels
        uicontrol(ph(which),'Style','text',...
                        'String','0',...
                        'Units','normalized',...
                        'Position',[0.1 0.05 0.08 0.2]);
        uicontrol(ph(which),'Style','text',...
                        'String','0.5',...
                        'Units','normalized',...
                        'Position',[0.475 0.05 0.05 0.2]);
        uicontrol(ph(which),'Style','text',...
                        'String','1',...
                        'Units','normalized',...
                        'Position',[0.85 0.05 0.05 0.2]);
        %button to zero slider
        uicontrol(ph(which),'Style','pushbutton',...
                        'String','1/2',...
                        'Units','normalized',...
                        'Position',[0.0 0.25 0.1 0.7],...
                        'Callback',{@pb_callback});
        %display value of slider
        sth(which)=uicontrol(ph(which),'Style','text',...
                        'String',num2str(fp),...
                        'Units','normalized',...
                        'Position',[0.9 0.25 0.1 0.7]);
    end%function MakeSlider

    function MakeEditBox(which)
        ph(which)=uipanel(figure1,'Units','normalized',...
                'Position',[0.87 0.05 0.1 0.1]);%left bottom width height
        set(ph(which),'Title','(# states)/2');
        %control
        uicontrol(ph(which),'Style','edit',...
            'String',num2str(nstates),...
            'Max',1,'Min',0,...
            'Units','normalized',...
            'Position',[0.05 0.05 0.9 0.9],...
            'Callback',{@edb_callback});
        %labels
%         uicontrol(ph(which),'Style','text',...
%             'String',title,...
%             'Units','normalized',...
%             'Position',[x y 1/6 1/2]);
    end%function MakeEditBox

    function MakePlot
        TestCorder(2*nstates,fp,'Parent',axes1);
    end%function MakePlot



end %function IntGui