function [ f,gr ] = OptFunGradL( x,s,fp,w,varargin )
%OPTFUNL function and gradient for fmincon on Laplace transform
%   x = parameters (off diagonal matrix elements)
%   f = function to be minimised (-SNR(t))
%   G = df/dx

[Wp,Wm]=Params2Mats(x);

ev=ones(1,length(Wp));
q=Wp-Wm;
W=Wm+fp*q;
p=EqProb(W);

Zinv=ones(length(Wp))-W;
Zinvs=s*eye(length(Wp))+Zinv;

a=q*(Zinvs\w);
c=(p*q)/Zinvs;

f=-2*fp*(1-fp)*p*a;

%dA(s)/dq_(i,j)
dAdq = ((Zinvs\w)*p)';
dAdq=dAdq-diag(dAdq)*ev;
%dA(s)/dW^F_(i,j)
dAdW = ((Zinv\a)*p)' + ((Zinvs\w)*c)';
dAdW=dAdW-diag(dAdW)*ev;
%dA/dWp_(i,i+1)+dA/dWm_(n-i+1,n-i)
dAdWp=dAdq+fp*dAdW;
dAdWm=-dAdq+(1-fp)*dAdW;

gr = -Mats2Params(dAdWp,dAdWm);


% f=-real(crq(Wp,Wm,fp,w));
% f=-SNRarea( Wp, Wm, fp, w );
end

