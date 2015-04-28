function [ qp,qm ] = MakeHomq( qv,fp )
%[qp_eff,qm_eff]=MAKEHOMQ(qv,fp) effective nearest neighbour transition rates with
%activity independent plasticity
%   qv=[qp qm qhp qhm]
%   qp_eff = qp + qhp/28fp
%   qm_eff = qm + qhm/28fp

n=length(qv)/4;

%activity dependent plasticity
qp=qv(1:n);
qm=qv(n+1:2*n);
%acticity independent plasticity
qhp=qv(2*n+1:3*n);
qhm=qv(3*n+1:end);

qp=qp + qhp/(2*fp);
qm=qm + qhm/(2*(1-fp));


end

