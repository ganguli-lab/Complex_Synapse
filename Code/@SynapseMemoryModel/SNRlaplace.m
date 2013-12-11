function [ A ] = SNRlaplace( obj,s )
%A=SynapseMemoryModel.SNRLAPLACE Laplace transform of SNR curve for complex synapse (cts time)
%   A(s) = int exp(-s*t)*SNR(t) dt
%   s  = parameter of Laplace transform
%   WP = potentiation transition rates
%   WM = depression transition rates
%   FP = Fraction of potentiation transitions
%   w  = Weights of states (+/-1)

if isscalar(s)
    Zinv=obj.GetZinv;
    A=(2*obj.fp*(1-obj.fp)) * sum( (Zinv\obj.GetEnc) * ((s*eye(length(obj.w))+Zinv)\obj.w));
else
    A=zeros(size(s));
    for i=1:numel(s)
        A(i)=SNRlaplace(s(i),Wp,Wm,fp,w);
    end
end


end

