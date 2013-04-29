function [ ih ] = ImTrans( Wp,Wm,Parent)
%IMTRANS Summary of this function goes here
%   Detailed explanation goes here

if ~exist('Parent','var')
    Parent=figure('WindowStyle','docked');
end
if ~exist('axp','var')
    axp=subplot(1,2,1,'Parent',Parent);
end
if ~exist('axm','var')
    axm=subplot(1,2,2,'Parent',Parent);
end

ih(1)=imagesc(Wp+eye(length(Wp)),'Parent',axp);
title(axp,'M^{pot}');
xlabel(axp,'To');
ylabel(axp,'From');
colorbar('peer',axp);

ih(2)=imagesc(Wm+eye(length(Wm)),'Parent',axm);
title(axm,'M^{dep}');
xlabel(axm,'To');
ylabel(axm,'From');
colorbar('peer',axm);


end

