function [ KTsgn,KTsat ] = NumLaplaceKT( chains,sym )
%NUMLAPLACEKT Summary of this function goes here
%   Detailed explanation goes here

KTsgn=zeros(size(chains));
KTsat=KTsgn;
nstates=length(chains(1).qv)/2 + 1;

for i=1:length(chains)
    
    if sym
        qp=chains(i).qv;
        qm=wrev(qp);
    else
        qp=chains(i).qv(1:nstates-1);
        qm=chains(i).qv(nstates:end);
    end
    [Wp,Wm]=MakeMultistate(qp,qm);
    
    Wp=Wp+eye(length(Wp));
    Wm=Wm+eye(length(Wp));
    
    KTsgn(i) = min( min(min(chains(i).KTp)), min(min(chains(i).KTm)) );
    
    KTsat(i) = max( max(max( abs(Wp.*chains(i).KTp) )), max(max( abs(Wm.*chains(i).KTm) )) );
end

end

