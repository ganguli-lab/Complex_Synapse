function [ Wp,Wm,w ] = MakeMultistate( qp,qm )
%[WP,WM,w]=MAKEMULTISTATE(QP,QM) Transition rates for Multistate toplogy
%   QP/QM=adjacent transiton rates for potentiation/depression, vector

error(CheckSize(qp,@isvector));
existsAndDefault('qm',wrev(qp));
error(CheckSize(qm,@(x) samesize(x,qp),'samesize(qp)'));
error(CheckSize(qp,@(x) mod(length(x),2)==1,'odd dim'));

Wp=StochastifyC(diag(qp,1));
Wm=StochastifyC(diag(qm,-1));

if nargout>2
w=[-ones((length(qp)+1)/2,1);ones((length(qp)+1)/2,1)];
end

end

