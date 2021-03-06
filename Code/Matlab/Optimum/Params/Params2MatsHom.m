function [ Wp,Wm,Q ] = Params2MatsHom( x )
%[Wp,Wm]=PARAMS2MATSHOM(x) convert vector of independent matrix elements to
%matrices
%   WP = potentiation transition rates
%   WM = depression transition rates
%   Q  = activity independent transition rates
%   x  = off diagonal elments, row major

error(CheckSize(x,@isvector));

n=length(x)/3;
n=round((1+sqrt(1+4*n))/2);

Wp=reshape(x(1:(n*(n-1))),n,n-1);
Wp=[reshape([zeros(1,n-1);Wp],1,[]) 0];
Wp=reshape(Wp,n,n);

Wm=reshape(x((n*(n-1)+1):2*n*(n-1)),n,n-1);
Wm=[reshape([zeros(1,n-1);Wm],1,[]) 0];
Wm=reshape(Wm,n,n);

Q=reshape(x((2*n*(n-1)+1):end),n,n-1);
Q=[reshape([zeros(1,n-1);Q],1,[]) 0];
Q=reshape(Q,n,n);

Wp=StochastifyC(Wp');
Wm=StochastifyC(Wm');
Q=StochastifyC(Q');

end

