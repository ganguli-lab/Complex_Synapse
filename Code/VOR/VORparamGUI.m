function VORparamGUI(varargin)
%VORparamGUI(...param,value,...) GUI for comparing learning curves for
%different parameters
%
%   builder_h = function handle that builds transition matrices (see folder
%               Builders, default: @SerialBuilder)
%   FontSize  = Fonst size for control labels (default: 16)

%initial params
persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='VORparamGUI';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addParameter('builder_h',@SerialBuilder);
    p.addParameter('FontSize',16);
%     p.addParameter('Normalise',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
end
p.parse(varargin{:});
r=p.Results;
% 
%check validity of params
assert(isa(r.builder_h,'function_handle'));

%Initialise data (put them in a global struct)
S.builder_h=r.builder_h;
S.numStates=10;
S.paramPot=0.5;
S.paramWT=0.5;
S.paramKO=0.5;
S.fpNorm=1/2;
S.fpInc=1/2;
S.fpDec=1/2;
S.t_inc=5;
S.t_dec=20;
%parameters to be controlled by sliders
slparams={'paramPot','paramWT','paramKO','fpNorm','fpInc','fpDec'};
%number of parameters to be controlled by edit-boxes
edparams={'numStates','t_inc','t_dec'};
%number of parameters to be controlled by sliders
numsl=length(slparams);
%number of parameters to be controlled by edit-boxes
numed=length(edparams);

%-------------------------------------------------------------------------
%Create figure
figure1 = figure('Units','normalized');

%Create axes
axes1 = axes('Parent',figure1,'OuterPosition',[0 0 0.5 1]);%left bottom width height

%where slider and edit-boxes will be
slpos={0.5,0.2,0.5,0.8};%left bottom width height
edpos={0.5,0,0.5,0.2};%left bottom width height

%-------------------------------------------------------------------------
%Make dummy handles
phs=zeros(1,numsl);
slh=phs;
sth=phs;
phe=zeros(1,numed);

%-------------------------------------------------------------------------
%Create param panels
for i=1:numsl;
    MakeSlider(i,slparams{i});
end
for i=1:numed;
    MakeEditBox(i,edparams{i})
end

%-------------------------------------------------------------------------
%Make initial plots
MakePlot;

%-------------------------------------------------------------------------

%  Callbacks for PertGui

    %callback for sliders
    function sl_callback(source,~,which,varname)
        S.(varname)=get(source,'Value');
        set(sth(which),'String',num2str(S.(varname)));
        MakePlot;
    end 

    %callback for push-buttons in slider panels
    function pb_callback(~,~,which,varname,val)
        set(slh(which),'Value',val);
        sl_callback(slh(which),0,which,varname);
    end

    %callback for edit-boxes in slider panels
    function edsl_callback(source,~,which,varname)
        val=str2double(get(source,'String'));
        set(slh(which),'Value',val);
        sl_callback(slh(which),0,which,varname);
    end

    %callback for edit-boxes
    function edb_callback(source,~,varname)
%         S.(varname)=floor(str2double(get(source,'String')));
        S.(varname)=str2double(get(source,'String'));
        MakePlot;
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

    function MakeSlider(which,varname)
        %slider values for min, max, middle (used for push-button)
        minval=0;
        maxval=1;
        midval=0.5;
        % midval=S.(varname);
        %make panel
        phs(which)=uipanel(figure1,'Units','normalized',...
                'Position',CalcPosVert(which,numsl,slpos{:}));%left bottom width height
        set(phs(which),'Title',varname,'FontSize',r.FontSize);
        %slider
        slh(which)=uicontrol(phs(which),'Style','slider',...
                        'Max',maxval,'Min',minval,'Value',S.(varname),...
                        'SliderStep',[0.002 0.1],...
                        'Units','normalized',...
                        'Position',[0.1 0.25 0.8 0.7],...
                        'Callback',{@sl_callback,which,varname});
        %slider labels
        uicontrol(phs(which),'Style','text',...
                        'String',num2str(minval),...
                        'Units','normalized',...
                        'FontSize',r.FontSize,...
                        'Position',[0.1 0.05 0.08 0.2]);
        uicontrol(phs(which),'Style','text',...
                        'String',num2str(midval),...
                        'Units','normalized',...
                        'FontSize',r.FontSize,...
                        'Position',[0.1+0.75*(midval-minval)/(maxval-minval) 0.05 0.05 0.2]);
        uicontrol(phs(which),'Style','text',...
                        'String',num2str(maxval),...
                        'Units','normalized',...
                        'FontSize',r.FontSize,...
                        'Position',[0.85 0.05 0.05 0.2]);
        %button to set slider to middle value
        uicontrol(phs(which),'Style','pushbutton',...
                        'String',num2str(midval),...
                        'Units','normalized',...
                        'FontSize',r.FontSize,...
                        'Position',[0.0 0.25 0.1 0.7],...
                        'Callback',{@pb_callback,which,varname,midval});
        %display value of slider
        sth(which)=uicontrol(phs(which),'Style','edit',...
                        'String',num2str(get(slh(which),'Value')),...
                        'Units','normalized',...
                        'FontSize',r.FontSize,...
                        'Position',[0.9 0.25 0.1 0.7],...
                        'Callback',{@edsl_callback,which,varname});
    end%function MakeSlider

    function MakeEditBox(which,varname)
        phe(which)=uipanel(figure1,'Units','normalized',...
                'Position',CalcPosHorz(which,numed,edpos{:}));%left bottom width height
        set(phe(which),'Title',varname,'FontSize',r.FontSize);
        %control
        uicontrol(phe(which),'Style','edit',...
            'String',num2str(S.(varname)),...
            'Max',1,'Min',0,...
            'Units','normalized',...
                        'FontSize',r.FontSize,...
            'Position',[0.05 0.05 0.9 0.9],...
            'Callback',{@edb_callback,varname});
        %labels
%         uicontrol(phe(which),'Style','text',...
%             'String',title,...
%             'Units','normalized',...
%             'Position',[x y 1/6 1/2]);
    end%function MakeEditBox

    function MakePlot
        vorexptobj=VORbuilder(S.builder_h,S.numStates,S.paramPot,S.paramWT,S.paramKO,S.fpNorm,S.fpInc,S.fpDec,S.t_inc,S.t_dec);
        vorexptobj.PlotLearnS('Parent',axes1,'LineWidth',2);
    end%function MakePlot



end %function IntGui