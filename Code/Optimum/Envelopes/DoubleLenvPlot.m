function [ ph,yl ] = DoubleLenvPlot( env,AenvSingle )
%[ph,yl]=DOUBLELENVPLOT(env) plot constrained memory curve envelope at times env.tau
%   env is struct made by NumLaplaceBndDouble
%   constraint is A(enc.sc)=Ac
%   


ph=loglog(env.tau,AenvSingle./env.tau,'g',env.tau,env.Aenv,'r',1/env.sc,env.sc*env.Ac,'rd');
hold on;
yl=ylim;

end

