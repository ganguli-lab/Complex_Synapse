function TwoModelPlot( chains,env,ind,varargin )
%TWOMODELPLOT Summary of this function goes here
%   Detailed explanation goes here


loglog(env.tau,[chains.s].*[chains.A],'g',env.tau,env.SNRbenv,'r',1/env.sc,env.sc*env.Ac,'rd',varargin{:});
hold on
yl=ylim;
if exist('ind','var') && ~isempty(ind)
    loglog(env.tau,env.chains(ind).snrb,'b',varargin{:});
    line([1 1]*env.tau(ind),yl,'Color','k',varargin{:});
    line([1 1]*env.tau(env.chains(ind).sinds(1)),yl,'Color','k','LineStyle',':',varargin{:});
    line([1 1]*env.tau(env.chains(ind).sinds(2)),yl,'Color','k','LineStyle',':',varargin{:});
end
hold off;

end

