function [ ih ] = ImTrans( Wp,Wm,varargin)
%ih=IMTRANS(Wp,Wm,Parent,axp,axm) image transition probs
%   ih = vector of image handles
%   Parent = figure handle

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='ImTrans';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addRequired('Wp',@(x)validateattributes(x,{'numeric'},{'2d','square'},'ImTrans','Wp',1));
    p.addRequired('Wm',@(x)validateattributes(x,{'numeric'},{'2d','square'},'ImTrans','Wm',2));
    p.addOptional('Parent',[],@(x)validateattributes(x,{'cell'},{},'ImTrans','Parent',3));
    p.addOptional('axp',[],@(x)validateattributes(x,{'cell'},{},'ImTrans','axp',4));
    p.addOptional('axm',[],@(x)validateattributes(x,{'cell'},{},'ImTrans','axm',5));
    p.addParameter('CLim',[0 1],@(x) validateattributes(x,{'numeric'},{'row','size',[1 2]},'ImTrans','CLim'));
end
p.parse(Wp,Wm,varargin{:});
r=p.Results;

if isempty(r.Parent) && (isempty(r.axp) || isempty(r.axm))
    r.Parent=figure('WindowStyle','docked');
end
if isempty(r.axp)
    r.axp=subplot(1,2,1,'Parent',r.Parent);
end
if isempty(r.axm)
    r.axm=subplot(1,2,2,'Parent',r.Parent);
end


cla(r.axp);
ih(1)=imagesc(Wp+eye(size(Wp)),'Parent',r.axp,p.Unmatched,r.CLim);
title(r.axp,'M^{pot}'); 
xlabel(r.axp,'To');
ylabel(r.axp,'From');
colorbar('peer',r.axp);

cla(r.axm);
ih(2)=imagesc(Wm+eye(size(Wm)),'Parent',r.axm,p.Unmatched,r.CLim);
title(r.axm,'M^{dep}');
xlabel(r.axm,'To');
ylabel(r.axm,'From');
colorbar('peer',r.axm);


end

