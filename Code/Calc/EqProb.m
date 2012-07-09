function [ pinf ] = EqProb( W )
%EQPROB(W) equlibrium distribution of Markov chain (cts time)
%   W = transition rates

error(CheckSize(W,@ismat));%matrix
error(CheckSize(W,@issquare));%square

pinf = ones(1,size(W,1))/(ones(size(W)) - W);

% [v,qb]=eig(-W');
% qb=diag(qb);
% [~,ix]=sort(qb);
% v=conj(v(:,ix))';
% pinf=v(1,:);
% pinf=pinf/sum(p);

end

