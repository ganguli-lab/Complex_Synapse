function [ tf,ix ] = istransient( W,thresh,varargin )
%[tf,ix]=ISTRANSIENT(W,thresh) Is Markov chain transient?
%   tf = are there any transient states?(true/false
%   ix = indices of transient states
%   W = transition rate matrix
%   thresh = threshold for transience, sum_i(i~=j) W_ij < thresh, (default: 1e-5)

error(CheckSize(W,@ismat));%matrix
error(CheckSize(W,@issquare));%square
existsAndDefault('thresh',1e-5);
error(CheckSize(thresh,@isscalar));%matrix

UseP=false;
varargin=assignApplicable(varargin);

if UseP
    p=EqProb(W);
    tf= p <thresh;
else
    tf = sum(W,1)-diag(W)' < thresh;
end

ix=find(tf);
tf=any(tf);


end

