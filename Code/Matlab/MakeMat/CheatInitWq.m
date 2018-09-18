function [ qc ] = CheatInitWq( W,q,w )
%qc=CHEATINIT(W,q,w) produce q-matrix that puts all initial SNR at ends
%   W = forgetting transition rates: f^+W^+ + f^-W^-
%   q = encoding transition rates: W^+ - W^-
%   w = Weights of states (+/-1)
%   qc = cheating version of q

error(CheckSize(W,@ismat));%matrix
error(CheckSize(W,@issquare));%square
error(CheckSize(q,@(x)samesize(W,x),'samesize(W)'));%also square matrix of same size
error(CheckSize(w,@iscol));


p=EqProb(W);
deta=DeltaEta(W,w);

[~,mini]=max(deta);
[~,maxi]=min(deta);

qc=zeros(size(q));

qc(:,mini)=-1;
qc(:,maxi)=1;

qc=qc*(p*q*w)/2;



end

