function [ Wp,Wm,w ] = CliqueChain( n )
%CLIQUECHAIN Summary of this function goes here
%   Detailed explanation goes here

assert(isscalar(n));
assert(mod(n,6)==0);

w=[-ones(n/2,1);ones(n/2,1)];

Wp=zeros(n);

nn=n/3;

Wp=AddFan(Wp,1:nn,nn+1,1);

Wp=AddChain(Wp,(nn+1):(nn+2),1);
Wp=AddChain(Wp,(nn+2):(2*nn-1),1);
Wp=AddChain(Wp,(2*nn-1):(2*nn),1);

Wp=AddFan(Wp,2*nn,(2*nn+1):n,1/nn);

Wp=AddClique(Wp,(2*nn+1):n,1/(nn-1));

M=max(-diag(Wp));
Wp=Wp/M;

Wm=rot90(Wp,2);


end

