function [ h,yl ] = PlotEnvs( t,n,varargin )
%h=PLOTENVS(t,n,...) Plot of SNR envelopes, proven and conjectured
%   t = time (row vector)
%   n = # states
%   h = plot handles
%   yl= y limits
%   ... Param/Value pairs
%       Format = do we format the graph? (default=true)

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='PlotEnvs';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addRequired('t',@(x)validateattributes(x,{'numeric'},{'row','nonnegative'},'PlotEnvs','t1',1));
    p.addRequired('n',@(x)validateattributes(x,{'numeric'},{'scalar','positive'},'PlotEnvs','n',2));
    p.addOptional('extraArgs',{},@(x)validateattributes(x,{'cell'},{},'PlotEnvs','extraArgs',3));
    p.addParameter('Format',true,@(x) validateattributes(x,{'logical'},{'scalar'},'PlotEnvs','Format'));
    p.addParameter('Parent',gca,@(x) validateattributes(x,{'numeric'},{'scalar'},'PlotEnvs','Parent'));
    p.addParameter('Constr3',true,@(x) validateattributes(x,{'logical'},{'scalar'},'PlotEnvs','Constr3'));
end
p.parse(t,n,varargin{:});
r=p.Results;



[env,bnds]=SNRenvelope(r.t,r.n);
[env2,bnds2]=SNRenvelope2(r.t,r.n);
if r.Constr3
    h=plot(r.t,env2,'g',r.t,env,'g--','LineWidth',3,'Parent',r.Parent,r.extraArgs{:});
else
    h=plot(r.t,env2,'g','LineWidth',3,'Parent',r.Parent,r.extraArgs{:});
end
hold(r.Parent,'on');
yl=[env(end) env2(1)];
if r.Constr3
    lh1=line([1;1]*bnds',yl'*ones(1,length(bnds)),'Color','k','LineWidth',1.5,'LineStyle',':','Parent',r.Parent,r.extraArgs{:});
else
    lh1=[];
end
lh2=line([1;1]*bnds2',yl'*ones(1,length(bnds2)),'Color','k','LineWidth',1.5,'LineStyle','--','Parent',r.Parent,r.extraArgs{:});
if r.Format
    set(r.Parent,'XScale','log','YScale','log','FontSize',16)
    xlim(r.Parent,[r.t(1) r.t(end)])
    ylim(r.Parent,yl)
    xlabel(r.Parent,'Time','FontSize',16)
    ylabel(r.Parent,'SNR','FontSize',16);
end
h=[h;lh1;lh2];
end

