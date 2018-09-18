function [ mh,lh ] = DoubleModelPlot( env,i,yl,varargin )
%[mh,lh]=DOUBLEMODELPLOT(env,i,...) plot memory curve for model in
%env.mats(i).modelobj at times env.tau
%   env is struct made by NumLaplaceBndDouble
mh=loglog(env.tau,env.mats(i).modelobj.SNRrunAve(env.tau),'b',varargin{:});
lh=line([1 1]/env.mats(i).s,yl,'Color','k',varargin{:});


end

