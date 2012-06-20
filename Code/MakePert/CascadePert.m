function [ Pertp,Pertm ] = CascadePert( W )
%[PERTP,PERTM] = SHORTCUTPERT(W,start,len) Perturbation that interpolates
%between multistate and cascade topologies
%   PERTP,PERTM = Perturbation of transition rates for potentiation,depression
%   W = matrix we're permuting
%   start,len = start,len of states we're shortcutting

error(CheckSize(W,@ismat));
error(CheckSize(W,@issquare));

n=length(W);
error(CheckValue(n,@(x) mod(x,2)==0));
nn=n/2;

%preallocate
Pertp=zeros(n);
Pertm=Pertp;

%calculate size of each perturbation. ensures multistate->0 at end
p=EqProb(W);
qp=[0;diag(W,1);0];%multistate trans rates
qp=qp.*[0;p'];%correct for prefactor in pert
qp=diff(qp);%size of pert needed
qm=[0;diag(W,-1);0];
qm=qm.*[p';0];
qm=-diff(qm);


for i=2:nn
    dpertp=ShortcutPert(W,nn-i+1,i);
    Pertp=Pertp+qp(nn-i+1)*dpertp;
    [~,dpertm]=ShortcutPert(W,nn,i);
    Pertm=Pertm+qm(nn+i)*dpertm;
end

end

