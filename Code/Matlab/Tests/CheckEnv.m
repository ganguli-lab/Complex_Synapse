function [ tf,S,t,env ] = CheckEnv( Wp,Wm,fp,w )
%[TF,S]=CHECKENV(WP,WM,FP,w) Does SNR exceed envelope
%   WP = potentiation transition rates
%   WM = depression transition rates
%   FP = Fraction of potentiation transitions
%   w = Weights of states (+/-1)
%   TF = does it exceed envelope? logical{0,1}
%   S = SNR curve

error(CheckSize(Wp,@ismat));%matrix
error(CheckSize(Wp,@issquare));%square
error(CheckSize(Wm,@(x)samesize(Wp,x),'samesize(Wp)'));%also square matrix of same size
error(CheckSize(fp,@isscalar));
error(CheckValue(Wp,@(x) inrange(fp,0,1),'inrange(0,1)'));%fp in [0,1]
error(CheckSize(w,@iscol));
error(CheckValue(w,@(x) all(x.^2==1),'all w = +/-1'));


n=length(Wp);
t=0:0.5:(2*(n-1)^2);

S=SNRcurve(t,Wp,Wm,fp,w);
% Su=SNRcurve(t,Wpu,Wmu,fp,w);
env=SNRenvelope(t,n);

tf=any(S(2:end)>env(2:end));


end

