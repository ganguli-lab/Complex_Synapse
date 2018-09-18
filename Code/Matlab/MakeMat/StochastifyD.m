function [ newW ] = StochastifyD( W )
%NEWW=STOCHASTIFY(W) turn into discrete time stochastic matrix
%   Scales rows so that sum_j W_ij = 1

validateattributes(W,{'numeric'},{'square'});
% assert(ismat(W));
% assert(issquare(W));

W(W<0 & W>-1e-7)=0;
newW=diag(1./sum(W,2))*W;


end

