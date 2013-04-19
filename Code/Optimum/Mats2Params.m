function [ x ] = Mats2Params( Wp,Wm )
%x=MATS2PARAMS(Wp,Wm) convert vector of independent matrix elements to
%matrices
%   WP = potentiation transition rates
%   WM = depression transition rates
%   x  = off diagonal elments, row major

error(CheckSize(Wp,@ismat));%matrix
error(CheckSize(Wp,@issquare));%square
error(CheckSize(Wm,@(x)samesize(Wp,x),'samesize(Wp)'));%also square matrix of same size

n=length(Wp);

x=reshape(Wp',[],1);
x(1:n+1:end)=[];

x=[x; reshape(Wm',[],1)];
x((n*(n-1)+1):(n+1):end)=[];


end

