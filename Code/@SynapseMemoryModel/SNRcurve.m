function [ S ] = SNRcurve( obj,t, varargin  )
%S=SynapseMemoryModel.SNRCURVE(T) SNR as function of time
%   T = time values

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='SynapseMemoryModel.SNRcurve';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addParameter('UseExpM',false,@(x)validateattributes(x,{'logical'},{'scalar'},'SynapseMemoryModel.SNRcurve','UseExpM'));
end
p.parse(varargin{:});
r=p.Results;


if r.UseExpM
    S=zeros(size(t));
    pr=obj.EqProb(p.Unmatched);
    q=obj.GetEnc;
    W=obj.GetWf;
    for i=1:numel(t)
        S(i) = pr*q*expm(W*t(i))*obj.w;
    end
    S=2*obj.fp*(1-obj.fp)*S;
else
    [ qa,ca ] = obj.Spectrum(p.Unmatched);
    S = (ca.*qa)'* exp(-qa*t);
end


end

