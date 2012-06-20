function CascadeGUI(varargin)

%initial params
nstates=4;
varargin=assignApplicable(varargin);
% 
%check validity of params
assert(isscalar(nstates));
assert(isint(nstates));



%Initialise data (put them in a global struct)
S.tmax=20;
S.tsteps=20;
S.neval=4;
S.fp=1/2;
S.pertval=0;
for i=1:(nstates/2)
    S.(['q' int2str(i)])=0.3;
end
%number of parameters to be controlled by sliders
numsl=2+nstates/2;
%number of parameters to be controlled by edit-boxes
numed=3;

%-------------------------------------------------------------------------
%Create figure
figure1 = figure('Units','normalized');

%Create axes
axes1 = axes('Parent',figure1,'OuterPosition',[0 0.5 0.5 0.5]);%left bottom width height
axes2 = axes('Parent',figure1,'OuterPosition',[0 0 0.5 0.5]);%left bottom width height

%where slider and edit-boxes will be
slpos={0.5,0.2,0.5,0.8};%left bottom width height
edpos={0.5,0.0,0.5,0.2};%left bottom width height

%-------------------------------------------------------------------------
%Make dummy handles
phs=zeros(1,numsl);
slh=phs;
sth=phs;
phe=zeros(1,numed);

%-------------------------------------------------------------------------
%Create param panels
MakeSlider(1,'fp');
for i=1:(nstates/2)
    MakeSlider(i+1,['q' int2str(i)]);
end
MakeSlider(numsl,'pertval');
MakeEditBox(1,'tmax')
MakeEditBox(2,'tsteps')
MakeEditBox(3,'neval')

%-------------------------------------------------------------------------
%Make initial plots
MakePlot;

%-------------------------------------------------------------------------

%  Callbacks for PertGui

    %callback for sliders
    function sl_callback(source,~,varname)
        S.(varname)=get(source,'Value');
        set(sth,'String',num2str(S.(varname)));
        MakePlot;
    end 

    %callback for push-buttons in slider panels
    function pb_callback(~,~,which,varname,val)
        set(slh(which),'Value',val);
        sl_callback(slh(which),0,varname);
    end

    %callback for edit-boxes
    function edb_callback(source,~,varname)
        S.(varname)=floor(str2double(get(source,'String')));
%         S.(varname)=str2double(get(source,'String'));
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
        set(phs(which),'Title',varname);
        %slider
        slh(which)=uicontrol(phs(which),'Style','slider',...
                        'Max',maxval,'Min',minval,'Value',S.(varname),...
                        'SliderStep',[0.02 0.2],...
                        'Units','normalized',...
                        'Position',[0.1 0.25 0.8 0.7],...
                        'Callback',{@sl_callback,varname});
        %slider labels
        uicontrol(phs(which),'Style','text',...
                        'String',num2str(minval),...
                        'Units','normalized',...
                        'Position',[0.1 0.05 0.08 0.2]);
        uicontrol(phs(which),'Style','text',...
                        'String',num2str(midval),...
                        'Units','normalized',...
                        'Position',[0.1+0.75*(midval-minval)/(maxval-minval) 0.05 0.05 0.2]);
        uicontrol(phs(which),'Style','text',...
                        'String',num2str(maxval),...
                        'Units','normalized',...
                        'Position',[0.85 0.05 0.05 0.2]);
        %button to set slider to middle value
        uicontrol(phs(which),'Style','pushbutton',...
                        'String',num2str(midval),...
                        'Units','normalized',...
                        'Position',[0.0 0.25 0.1 0.7],...
                        'Callback',{@pb_callback,which,varname,midval});
        %display value of slider
        sth(which)=uicontrol(phs(which),'Style','text',...
                        'String',num2str(get(slh(which),'Value')),...
                        'Units','normalized',...
                        'Position',[0.9 0.25 0.1 0.7]);
    end%function MakeSlider

    function MakeEditBox(which,varname)
        phe(which)=uipanel(figure1,'Units','normalized',...
                'Position',CalcPosHorz(which,numed,edpos{:}));%left bottom width height
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
        %make Markov chains
        q=zeros(1,nstates-1);
        for j=1:(nstates/2)
            q(j)=S.(['q' int2str(j)]);
            q(nstates-j)=q(j);
        end
        
        [Wp,Wm,w]=CascadeMSinterp(q,q,S.fp,S.pertval);
        
        %calc decay rates and area of eigenmodes
        [qa,ca]=Spectrum(Wp,Wm,S.fp,w);
        qa(1)=[];
        ca(1)=[];
        
        %calc SNR curve
        t=0:(S.tmax/S.tsteps):S.tmax;
        snr=SNRcurve(t,Wp,Wm,S.fp,w);
        
        %calc eigenmodes
        [u,d]=eig(S.fp*Wp+(1-S.fp)*Wm);
        [~,ix]=sort(-diag(d));
        u=u(:,ix);
        
        %plot SNR curve and time-constant/initial-SNR of eigenmodes
        plot(axes1,t,real(snr));
        hold(axes1,'all');
        plot(axes1,real(1./qa),real(ca.*qa),'d');
        xlabel(axes1,'Time');
        ylabel(axes1,'SNR');
        hold(axes1,'off');
        
        %plot eigenmodes
        plot(axes2,real(u(:,1:S.neval)));
        legend(axes2,int2str((1:S.neval)'));
        xlabel(axes2,'State');
        ylabel(axes2,'Eigenfunction');
    end%function MakePlot

end %function IntGui