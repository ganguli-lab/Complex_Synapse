function MakeVORfigs( prefix,paramWT,paramKn,df,T_train,T_pre,PrintView,varargin )
%MAKEVORFIGS(prefix,paramWT,paramKn,df,T_train,T_pre,PrintView) comparison of VOR learning
%curves for a particular synaptic model
%   paramWT   = parameter used for wild-type (pot&dep) and knockout (pot)
%   paramKn   = parameter used for knockout (dep)
%   df        = change in f+ vor gain increase learning
%   n         = number of internal synaptic states
%   modelname = name of synaptic model (multistate or cascade)


DoPrint= strcmpi(PrintView,'Print');

numpts=100;
close('all');

if DoPrint
    fig=figure('PaperPositionMode','auto','Position',[60 60 1000 1000]);
    figs=figure('PaperPositionMode','auto','Position',[60 60 1000 1000]);
    altfig(1)=figure('PaperPositionMode','auto','Position',[60 60 1000 500]);
    altfig(2)=figure('PaperPositionMode','auto','Position',[60 60 1000 500]);
    fig3(1)=figure('PaperPositionMode','auto','Position',[60 60 600 500]);
    fig3(2)=figure('PaperPositionMode','auto','Position',[60 60 600 500]);
    fig3(3)=figure('PaperPositionMode','auto','Position',[60 60 600 500]);
    fig3(4)=figure('PaperPositionMode','auto','Position',[60 60 600 500]);
else
    fig=figure('WindowStyle','docked','PaperPositionMode','auto');
    figs=figure('WindowStyle','docked','PaperPositionMode','auto');
    altfig=figure('WindowStyle','docked','PaperPositionMode','auto');
    fig3=figure('WindowStyle','docked','PaperPositionMode','auto');
end

if strncmpi(prefix,'cascade',3)
    modelname='cascade';
    n=10;
elseif strncmpi(prefix,'multistate',3)
    modelname='multistate';
    n=10;
elseif strncmpi(prefix,'binary',3)
    modelname='multistate';
    n=2;
elseif strncmpi(prefix,'pooled',3)
    modelname='pooled';
    n=10;
end



varargin=assignApplicable(varargin);

t=0:(T_train+T_pre)/numpts:(T_train+T_pre);
figure(fig);

MakeVORcomparison(paramWT,paramKn,df,n,modelname,'t',t,'tchange',T_pre,...
    'Superpose',false,'fig',fig,'altfig',altfig,'fig3',fig3,varargin{:})

t=0:T_train/numpts:T_train;
t=[t T_pre+t(t>T_train-T_pre)];

MakeVORcomparison(paramWT,paramKn,df,n,modelname,'t',t,'tchange',T_pre,...
    'Superpose',true,'fig',figs,varargin{:})

if DoPrint
    print(fig,[prefix '_learn.eps'],'-depsc');
    print(figs,[prefix '_learnS.eps'],'-depsc');
    print(altfig(1),[prefix '_eq_WT.eps'],'-depsc');
    print(altfig(2),[prefix '_eq_KO.eps'],'-depsc');
    print(fig3(1),[prefix '_pr_WT_nopre.eps'],'-depsc');
    print(fig3(2),[prefix '_pr_WT_pre.eps'],'-depsc');
    print(fig3(3),[prefix '_pr_KO_nopre.eps'],'-depsc');
    print(fig3(4),[prefix '_pr_KO_pre.eps'],'-depsc');
    close('all');
end


end

