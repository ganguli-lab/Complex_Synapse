function [ pr ] = MSeqProb( qp,qm,fp )
%pr=MSEQPROB(qp,qm,fp) equilibrium distribution for Multistate/Serial topology
%   pr = equilibrium distribution
%   qp = nearest neighbour transition rate for pot
%   qm = nearest neighbour transition rate for dep
%   fp = fraction of events that are pot

pr=(fp*qp)./((1-fp)*qm);
pr=[1 cumprod(pr)];
pr=pr/sum(pr);


end

