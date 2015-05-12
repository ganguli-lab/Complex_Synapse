function [ f,gr ] = OptFunGradHomL( x,s,fp,w,varargin )
%OPTFUNGRADHOML function and gradient for fmincon on Laplace transform
%   x = parameters (off diagonal matrix elements)
%   f = function to be minimised (-SNR(t))
%   G = df/dx

[Wp,Wm,Q]=Params2MatsHom(x);

Wp=Wp + 1/(2*fp) * Q;
Wm=Wm + 1/(2*(1-fp)) * Q;

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
dAdQ=dAdWp/(2*fp) + dAdWm/(2*(1-fp));

gr = -Mats2ParamsHom(dAdWp,dAdWm,dAdQ);


% f=-real(crq(Wp,Wm,fp,w));
% f=-SNRarea( Wp, Wm, fp, w );
end

