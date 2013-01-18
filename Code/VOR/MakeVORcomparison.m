function MakeVORcomparison( paramWT,paramKn,df,n,modelname,varargin )
%MAKEVORCOMPARISON(paramp,paramm,df,n,modelname) comparison of VOR learning
%curves for a particular synaptic model
%   paramWT   = parameter used for wild-type (pot&dep) and knockout (pot)
%   paramKn   = parameter used for knockout (dep)
%   df        = change in f+ vor gain increase learning
%   n         = number of internal synaptic states
%   modelname = name of synaptic model (multistate or cascade)


t=0:1:80;
tchange=40;
LinWeights=false;
gradq=0;
altfig=[];
fig=gcf;

varargin=assignApplicable(varargin);

if n>2
    dq=0:1/(n-2):1;
else
    dq=zeros(1,n-1);
end

if strcmpi(modelname,'multistate')
    [Wp,Wm,w]=MakeSMS(paramWT*(ones(1,n-1)+gradq*dq));
    [~,WmKn]=MakeSMS(paramKn*(ones(1,n-1)+gradq*dq));
elseif strcmpi(modelname,'cascade')
    [Wp,Wm,w]=CascadeOriginal(paramWT,paramWT,n);
    [~,WmKn]=CascadeOriginal(paramWT,paramKn,n);
else
    error(['invalid model name: ' modelname]);
end

if LinWeights
    w=(-1:2/(n-1):1)';
end

clf(fig);
set(fig,'WindowStyle','docked')
Parent=axes('Parent',fig);
if isempty(altfig)
    altfig=figure('WindowStyle','docked');
end

VORcomparison(Wp,Wm,w,t,df,tchange,Wp,WmKn,'LineWidth',2,'Parent',Parent,'altfig',altfig,varargin{:});

end

