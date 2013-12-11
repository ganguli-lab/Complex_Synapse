function [ area ] = SNRarea( obj )
%AREA=SynapseMemoryModel.SNRAREA Area under SNR curve for complex synapse (cts time)
%   WP = potentiation transition rates
%   WM = depression transition rates
%   FP = Fraction of potentiation transitions
%   w  = Weights of states (+/-1)



% n=size(Wp,1);
% assert(mod(n,2)==0)
% 
% q=Wp-Wm;
% w=ones(n,1);
% w(1:(n/2))=-1;

Zinv=obj.GetZinv;

area = 2*obj.fp *(1-obj.fp) * sum( (Zinv \ obj.GetEnc) * (Zinv \ obj.w) );
% area =  sum( (Zinv \ q) * (Zinv \ w) );


end

