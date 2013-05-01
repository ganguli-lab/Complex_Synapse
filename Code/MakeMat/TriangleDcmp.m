function [ Wp,Wm,neww ] = TriangleDcmp( W,fp,w,t )
%[WP,WM,neww]=TRIANGLEDCMP(W,FP,w) Extract W^+,W^- from W, consistent with
%upper/lower triangularity assocated with w.
%   Detailed explanation goes here
%   W = transition rates
%   FP = Fraction of potentiation transitions
%   w = Weights of states (+/-1), default=[-1;...;-1;1;...;1].
%   WP = potentiation transition rates
%   WM = depression transition rates
%   neww = w after sorting
%

error(CheckSize(W,@ismat));%matrix
error(CheckSize(W,@issquare));%square
error(CheckSize(fp,@isscalar));
error(CheckValue(fp,@(x) inrange(x,0,1),'inrange(0,1)'));%fp in [0,1]


if exist('w','var')
    %check w is valid
    assert(iscol(w));%row
    assert(length(w)==length(W));%same size
    assert(all(abs(w)==1));%+/-1
    %
    if exist('t','var')
        [W,neww]=SortByWt(W,w,t);
    else
        [W,neww]=SortByEta(W,w);
    end
else
    neww=[-ones(length(W)/2,1);ones(length(W)/2,1)];
end

Wp=triu(W)/fp;
Wm=tril(W)/(1-fp);
Wp=StochastifyC(Wp);
Wm=StochastifyC(Wm);


end

