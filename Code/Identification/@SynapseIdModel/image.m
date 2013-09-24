function imh=image(obj,axMInitial,axMM,varargin)
%imh=IMAGE(obj,axMInitial,axMM,varargin) image plots of probs
%   imh = image handles [initial,M{1},...]
%   obj = SynapseIdModel
%   axMInitial = axMes handle for initial (empty to not plot)
%   axMM       = [axMes handles] for M{:} (empty to not plot)

if length(obj.M)==2
    Mlabels={'pot','dep'};
else
    Mlabels=cellstr(int2str((1:length(obj.M))'))';
end
varargin=assignApplicable(varargin);

imh=[];
if ~isempty(axMInitial)
    imh=imagesc(obj.initial,'Parent',axMInitial);
    title(axMInitial,'Initial',varargin{:});
    xlabel(axMInitial,'State',varargin{:});
    set(axMInitial,'CLim',[0 1]);
end
for i=1:length(axMM)
    imh(end+1)=imagesc(obj.M{i},'Parent',axM(i));
    title(axM(i),['M^{' Mlabels{i} '}'],varargin{:});
    xlabel(axM(i),'To',varargin{:});
    ylabel(axM(i),'From',varargin{:});
    set(axM(i),'CLim',[0 1]);
end


end

