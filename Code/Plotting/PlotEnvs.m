function [ h,yl ] = PlotEnvs( t,n,varargin )
%h=PLOTENVS(t,n) Plot of SNR envelopes, proven and conjectured
%   Detailed explanation goes here

error(CheckSize(t,@isrow))
error(CheckSize(n,@isscalar))

format=true;
varargin=assignApplicable(varargin);

[env,bnds]=SNRenvelope(t,n);
[env2,bnds2]=SNRenvelope2(t,n);
h=plot(t,env2,'g',t,env,'g--','LineWidth',3,varargin{:});
hold on
yl=[env(end) env2(1)];
line([1;1]*bnds',yl'*ones(1,length(bnds)),'Color','k','LineWidth',1.5,'LineStyle',':',varargin{:});
line([1;1]*bnds2',yl'*ones(1,length(bnds2)),'Color','k','LineWidth',1.5,'LineStyle','--',varargin{:});
if format
    set(gca,'XScale','log','YScale','log')
    xlim([t(1) t(end)])
    ylim(yl)
    embiggen
    xlabel('Time')
    ylabel('SNR');
end

end

