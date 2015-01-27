function [ As,dAp,dAm ] = SNRlaplaceGrad( obj,s )
%[As,dAp,dAm]=obj.SNRLAPLACEGRAD(s) derivative of Laplace transform of
%SNR curve wrt matrix elemements
%   s       = value of Laplace transform parameter at which we evaluate
%   As      = Laplace transform of SNR curve at s
%   dAp(ij) = gradient of As wrt Wp(ij), i~=j
%   dAm(ij) = gradient of As wrt Wm(ij), i~=j
%   Assumes Wpm(ii) = -sum_(j~=i) Wpm(ij) enforced during all variations

ev=ones(1,obj.NumStates);
q=obj.GetEnc;
p=obj.EqProb;

Zinv=GetZinv(obj);
Zinvs=s*eye(obj.NumStates)+Zinv;

%this a_i is 2a_i/p_i from notes 
a=q*(Zinvs\obj.w);
%this c_k is c_k/p_k from notes 
c=(p*q)/Zinvs;

As=2*obj.fp*(1-obj.fp)*p*a;

%dA(s)/dq_(i,j)
dAdq = ((Zinvs\obj.w)*p)';
dAdq=dAdq-diag(dAdq)*ev;
%dA(s)/dW^F_(i,j)
dAdW = ((Zinv\a)*p)' + ((Zinvs\obj.w)*c)';
dAdW=dAdW-diag(dAdW)*ev;
%dA/dWp_(i,i+1)+dA/dWm_(n-i+1,n-i)
dAp=dAdq+obj.fp*dAdW;
dAm=-dAdq+(1-obj.fp)*dAdW;



end

