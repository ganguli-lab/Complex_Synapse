function DoubleLenvPlot( env,AenvSingle,ind,varargin )
%DOUBLELENVPLOT(env,AenvSingle,ind) plot constrained memory curve envelope at times env.tau
%   env is struct made by NumLaplaceBndDouble
%   constraint is A(enc.sc)=Ac
%   


loglog(env.tau,AenvSingle./env.tau,'g',env.tau,env.Aenv,'r',1/env.sc,env.sc*env.Ac,'rd',varargin{:});

hold on;

if exist('ind','var') && ~isempty(ind)
    loglog(env.tau,env.mats(ind).modelobj.SNRrunAve(env.tau),'b',varargin{:});
    yl=ylim;
    line([1 1]/env.mats(ind).s,yl,'Color','k',varargin{:});
end

hold off

end

