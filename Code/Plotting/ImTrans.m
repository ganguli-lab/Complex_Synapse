function [ ih ] = ImTrans( Wp,Wm,Parent,axp,axm)
%ih=IMTRANS(Wp,Wm,Parent,axp,axm) image transition probs
%   ih = vector of image handles
%   Parent = figure handle

if ~exist('Parent','var')
    Parent=figure('WindowStyle','docked');
end
if ~exist('axp','var')
    axp=subplot(1,2,1,'Parent',Parent);
end
if ~exist('axm','var')
    axm=subplot(1,2,2,'Parent',Parent);
end

cla(axp);
ih(1)=imagesc(Wp+eye(length(Wp)),'Parent',axp,[0 1]);
title(axp,'M^{pot}');
xlabel(axp,'To');
ylabel(axp,'From');
colorbar('peer',axp);

cla(axm);
ih(2)=imagesc(Wm+eye(length(Wm)),'Parent',axm,[0 1]);
title(axm,'M^{dep}');
xlabel(axm,'To');
ylabel(axm,'From');
colorbar('peer',axm);


end

