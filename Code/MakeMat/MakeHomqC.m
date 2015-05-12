function [ qp,qm ] = MakeHomqC( qv,fp )
%[qp_eff,qm_eff]=MAKEHOMQ(qv,fp) effective nearest neighbour transition rates with
%activity independent plasticity, for forgetting, not encoding.
%   qv=[qp qm qhp qhm]
%   qp_eff = qp + qhp/28fp
%   qm_eff = qm + qhm/28fp

n=length(qv)/4;

%activity dependent plasticity
qpi=qv(1:n);
qpd=qv(n+1:2*n);
qmi=qv(2*n+1:3*n);
qmd=qv(3*n+1:end);

qp=qpi + qmi*(1-fp)/fp;
qm=qmd + qpd*fp/(1-fp);


end

