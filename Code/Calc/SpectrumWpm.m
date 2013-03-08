function [ qa,ca ] = SpectrumWpm( Wp, Wm, fp, w )
%[QA,CA]=SPECTRUMWPM(WP,WM,FP,w) Weights and eigenvalues of SNR curve for complex synapse (cts time)
%   WP = potentiation transition rates
%   WM = depression transition rates
%   FP = Fraction of potentiation transitions
%   w = Weights of states (+/-1)
%   QA = eigenvalues (decay rate)
%   CA = weights (contribution to area)

error(CheckSize(Wp,@ismat));%matrix
error(CheckSize(Wp,@issquare));%square
error(CheckSize(Wm,@(x)samesize(Wp,x),'samesize(Wp)'));%also square matrix of same size
error(CheckSize(fp,@isscalar));
error(CheckValue(fp,@(x) inrange(fp,0,1),'inrange(0,1)'));%fp in [0,1]
error(CheckSize(w,@iscol));
error(CheckValue(w,@(x) all(x.^2==1),'all w = +/-1'));

n=size(Wp,1);
assert(mod(n,2)==0)

% w=ones(n,1);
% w(1:(n/2))=-1;

q=Wp-Wm;
W=fp*q + Wm;
[qa,ca]=SpectrumWq(W,q,fp,w);
% Zinv=ones(size(Wp)) - W;
% 
% [u,qa]=eig(-W);
% 
% % if rcond(u)<0.0001
% %     disp('bad u');
% % end
% % if rcond(Zinv)<0.0001
% %     disp('bad Z');
% % end
% ca = 2*fp*(1-fp) * (u\w) * sum((Zinv\q) * (Zinv\u), 1);
% 
% 
% qa=diag(qa);
% ca=diag(ca);
% 
% [qa,ix]=sort(qa);
% ca=ca(ix);
end

