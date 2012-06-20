function [ Wp,Wm,w ] = AllMtoP( n )
%ALLMTOP Summary of this function goes here
%   Detailed explanation goes here

error(CheckSize(n,@isscalar));%scalar
error(CheckValue(n,@(x) mod(x,2)==0));%even

nn=n/2;
Wp=Stochastify([zeros(nn), ones(nn); zeros(nn), zeros(nn)]/nn);

if nargout>1
Wm=rot90(Wp,2);
end

if nargout>2
w=[-ones(nn,1);ones(nn,1)];
end


end

