function [ S ] = SNRcurve( obj,t, varargin  )
%S=SynapseMemoryModel.SNRCURVE(T) SNR as function of time
%   T = time values


% S = gmdmp(ca.*qa, 1, exp(-outer(qa,t,true)), 1);
UseExpM=false;
varargin=assignApplicable(varargin);

if UseExpM
    S=zeros(size(t));
    p=obj.EqProb;
    q=obj.GetWf;
    W=obj.GetWf;
    for i=1:numel(t)
        S(i) = p*q*expm(W*t(i))*obj.w;
    end
    S=2*fp*(1-fp)*S;
else
    [ qa,ca ] = obj.Spectrum(varargin{:});
    S = (ca.*qa)'* exp(-qa*t);
end


end

