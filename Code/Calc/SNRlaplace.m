function [ A ] = SNRlaplace( s,Wp,Wm,fp,w )
%A=SNRLAPLACE(s,Wp,Wm,fp,w) Laplace transform of SNR curve for complex synapse (cts time)
%   A(s) = int exp(-s*t)*SNR(t) dt
%   s  = parameter of Laplace transform
%   WP = potentiation transition rates
%   WM = depression transition rates
%   FP = Fraction of potentiation transitions
%   w  = Weights of states (+/-1)
error(CheckSize(Wp,@ismat));%matrix
error(CheckSize(Wp,@issquare));%square
error(CheckSize(Wm,@samesize,'samesize(Wp)',Wp));%also square matrix of same size
error(CheckSize(fp,@isscalar));
error(CheckValue(fp,@inrange,'inrange(0,1)',0,1));%fp in [0,1]
error(CheckSize(w,@iscol));
error(CheckValue(w,@(x) all(x.^2==1),'all w = +/-1'));
error(CheckSize(w,@samelength,'samelength(Wp)',Wp));%same size

if isscalar(s)
    q=Wp-Wm;
    Zinv=ones(length(Wp))-Wm-fp*q;
    A=(2*fp*(1-fp)) * sum( (Zinv\q) * ((s*eye(length(Wp))+Zinv)\w));
else
    A=zeros(size(s));
    for i=1:numel(s)
        A(i)=SNRlaplace(s(i),Wp,Wm,fp,w);
    end
end


end

