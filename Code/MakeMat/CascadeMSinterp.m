function [Wp,Wm,w]=CascadeMSinterp(qp,qm,fp,lambda)
%[WP,WM,w]=CASCADEMSINTERP(QP,QM,FP,lambda) Interpolation between Multistate
%and Cascade topologies
%   QP,QM  = adjacent transiton rates for potentiation, row vectors
%   FP = Fraction of potentiation transitions, in [0,1]
%   lambda = how far to interpolate, in [0,1]
%   WP = potentiation transition rates
%   WM = depression transition rates
%   w  = Weights of states (+/-1)

error(CheckSize(qp,@isrow));
error(CheckSize(qp,@(x) mod(length(x),2)==1));
error(CheckSize(qm,@isrow));
error(CheckSize(qm,@(x)samesize(qp,x)));
error(CheckSize(fp,@isscalar));
error(CheckValue(fp,@(x) inrange(x,0,1)));
error(CheckSize(lambda,@isscalar));
error(CheckValue(lambda,@(x) inrange(x,0,1)));

n=length(qp)+1;
nn=n/2;


Wp=zeros(n);
Wm=Wp;

for i=1:nn-1
    Wp(i,i+1)=qp(i)*(1-lambda);
    Wp(nn+i,nn+i+1)=qp(nn+i);
    Wm(nn+i+1,nn+i)=qm(nn+i)*(1-lambda);
    Wm(i+1,i)=qm(i);
    Wp(i+1,nn+1)=(qp(i+1)-(1-fp)/fp*qm(i))*lambda;
    Wm(nn+i,nn)=(qm(nn+i-1)-fp/(1-fp)*qp(nn+i))*lambda;
end%for i

Wp(nn,nn+1)=qp(nn)-(1-fp)/fp*qm(nn-1)*lambda;
Wm(nn+1,nn)=qm(nn)-fp/(1-fp)*qp(nn+1)*lambda;
Wp(1,nn+1)=qp(1)*lambda;
Wm(n,nn)=qm(n-1)*lambda;

Wp=StochastifyC(Wp);
Wm=StochastifyC(Wm);

w=[-ones(nn,1);ones(nn,1)];

end%function


