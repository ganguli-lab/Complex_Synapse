function [ newW,neww,ix ] = SortByWt( W,w,t )
%[newW,neww,ix]=SORTBYWT(W,w) Put states in order of increasing exp(Wt)w
%   W = transition rates
%   w = Weights of states (+/-1)
%   t = time
%   ix=sort order
%   newW=W(ix,ix)
%   neww=w(ix)

error(CheckSize(W,@ismat));%matrix
error(CheckSize(W,@issquare));%square
error(CheckSize(w,@iscol));%row
error(CheckSize(w,@(x)length(x)==length(W),'same length as W'));%same size

deta=expm(W*t)*w;

[~,ix]=sort(deta,'ascend');

newW=W(ix,ix);
neww=w(ix);


end

