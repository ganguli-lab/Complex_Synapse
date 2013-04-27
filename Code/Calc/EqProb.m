function [ pinf ] = EqProb( W,varargin )
%EQPROB(W) equlibrium distribution of Markov chain (cts time)
%   W = transition rates

error(CheckSize(W,@ismat));%matrix
error(CheckSize(W,@issquare));%square

RCondThresh=1e-3;
varargin=assignApplicable(varargin);

Zinv=ones(length(W)) - W;

if rcond(Zinv)>RCondThresh
    pinf = ones(1,size(W,1))/(ones(size(W)) - W);
else
    [v,qb]=eig(-W');
    qb=diag(qb);
    [~,ix]=min(qb);
    pinf=v(:,ix).';
    pinf=pinf/sum(pinf);
end

end

