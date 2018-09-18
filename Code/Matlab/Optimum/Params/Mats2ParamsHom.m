function [ x ] = Mats2ParamsHom( Wp,Wm,Q )
%x=MATS2PARAMSHOM(Wp,Wm,Q) convert vector of independent matrix elements to
%matrices
%   WP = potentiation transition rates
%   WM = depression transition rates
%   Q  = activity independent transition rates
%   x  = off diagonal elments, row major

error(CheckSize(Wp,@ismat));%matrix
error(CheckSize(Wp,@issquare));%square
error(CheckSize(Wm,@(x)samesize(Wp,x),'samesize(Wp)'));%also square matrix of same size
error(CheckSize(Q,@issquare));%square
error(CheckSize(Q,@(x)samesize(Wp,x),'samesize(Wp)'));%also square matrix of same size

n=length(Wp);

x=reshape(Wp',[],1);
x(1:n+1:end)=[];

x=[x; reshape(Wm',[],1)];
x((n*(n-1)+1):(n+1):end)=[];

x=[x; reshape(Q',[],1)];
x((2*n*(n-1)+1):(n+1):end)=[];


end

