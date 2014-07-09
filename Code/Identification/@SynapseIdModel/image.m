function imh=image(obj,axMInitial,axM,varargin)
%imh=IMAGE(obj,axMInitial,axM,varargin) image plots of probs
%   imh = image handles [Initial,M{1},...]
%   obj = SynapseIdModel
%   axMInitial = axMes handle for Initial (empty to not plot)
%   axM        = [axMes handles] for M{:} (empty to not plot)


persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='SynapseIdModel.image';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addOptional('extraArgs',SynapseOptimset,@(x)validateattributes(x,{'cell'},{},'SynapseIdModel.image','extraArgs',4));
    p.addParameter('Mlabels',{'pot','dep'},@(x)validateattributes(x,{'cell'},{},'SynapseIdModel.image','Mlabels'));
%     p.addParameter('Normalise',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
end
p.parse(varargin{:});
r=p.Results;

if length(obj.M)~=length(r.Mlabels)
    r.Mlabels=cellstr(int2str((1:length(obj.M))'))';
end

imh=zeros(1,length(axM));
if ~isempty(axMInitial)
    imh=imagesc(obj.Initial,'Parent',axMInitial);
    title(axMInitial,'Initial',r.extraArgs{:});
    xlabel(axMInitial,'State',r.extraArgs{:});
    set(axMInitial,'CLim',[0 1]);
end
for i=1:length(axM)
    imh(i)=imagesc(obj.M{i},'Parent',axM(i));
    title(axM(i),['M^{' r.Mlabels{i} '}'],r.extraArgs{:});
    xlabel(axM(i),'To',r.extraArgs{:});
    ylabel(axM(i),'From',r.extraArgs{:});
    set(axM(i),'CLim',[0 1]);
end


end

