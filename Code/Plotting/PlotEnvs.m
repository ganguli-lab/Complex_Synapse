function [ h,yl ] = PlotEnvs( t,n,varargin )
%h=PLOTENVS(t,n,...) Plot of SNR envelopes, proven and conjectured
%   t = time (row vector)
%   n = # states
%   h = plot handles
%   yl= y limits
%   ... Param/Value pairs
%       Format = do we format the graph? (default=true)

error(CheckSize(t,@isrow))
error(CheckSize(n,@isscalar))

Format=true;
Parent=gca;
Constr3=true;
varargin=assignApplicable(varargin);

[env,bnds]=SNRenvelope(t,n);
[env2,bnds2]=SNRenvelope2(t,n);
if Constr3
    h=plot(t,env2,'g',t,env,'g--','LineWidth',3,'Parent',Parent,varargin{:});
else
    h=plot(t,env2,'g','LineWidth',3,'Parent',Parent,varargin{:});
end
hold(Parent,'on');
yl=[env(end) env2(1)];
if Constr3
    lh1=line([1;1]*bnds',yl'*ones(1,length(bnds)),'Color','k','LineWidth',1.5,'LineStyle',':','Parent',Parent,varargin{:});
else
    lh1=[];
end
lh2=line([1;1]*bnds2',yl'*ones(1,length(bnds2)),'Color','k','LineWidth',1.5,'LineStyle','--','Parent',Parent,varargin{:});
if Format
    set(Parent,'XScale','log','YScale','log','FontSize',16)
    xlim(Parent,[t(1) t(end)])
    ylim(Parent,yl)
    xlabel(Parent,'Time','FontSize',16)
    ylabel(Parent,'SNR','FontSize',16);
end
h=[h;lh1;lh2];
end

