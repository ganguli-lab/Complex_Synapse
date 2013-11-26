function [ ih ] = ImTrans( Wp,Wm,Parent,axp,axm,varargin)
%ih=IMTRANS(Wp,Wm,Parent,axp,axm) image transition probs
%   ih = vector of image handles
%   Parent = figure handle

if ~exist('Parent','var') || isempty(Parent)
    Parent=figure('WindowStyle','docked');
end
if ~exist('axp','var') || isempty(axp)
    axp=subplot(1,2,1,'Parent',Parent);
end
if ~exist('axm','var') || isempty(axm)
    axm=subplot(1,2,2,'Parent',Parent);
end

CLim=[0 1];
varargin=assignApplicable(varargin);

if ~isempty(CLim)
    varargin{end+1}=CLim;
end

cla(axp);
ih(1)=imagesc(Wp+eye(size(Wp)),'Parent',axp,varargin{:});
title(axp,'M^{pot}'); 
xlabel(axp,'To');
ylabel(axp,'From');
colorbar('peer',axp);

cla(axm);
ih(2)=imagesc(Wm+eye(size(Wm)),'Parent',axm,varargin{:});
title(axm,'M^{dep}');
xlabel(axm,'To');
ylabel(axm,'From');
colorbar('peer',axm);


end

