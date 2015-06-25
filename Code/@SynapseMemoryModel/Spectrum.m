function [ qa,ca ] = Spectrum( obj,varargin )
%[QA,CA]=SynapseMemoryModel.SPECTRUMWPM Weights and eigenvalues of SNR curve for complex synapse (cts time)
%   WP = potentiation transition rates
%   WM = depression transition rates
%   FP = Fraction of potentiation transitions
%   w = Weights of states (+/-1)
%   QA = eigenvalues (decay rate)
%   CA = weights (contribution to area)

% persistent p
% if isempty(p)
%     p=inputParser;
%     p.FunctionName='SynapseMemoryModel.Spectrum';
%     p.StructExpand=true;
%     p.KeepUnmatched=false;
%     p.addParameter('DegThresh',1e-3,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'},'SynapseMemoryModel.Spectrum','DegThresh'));
%     p.addParameter('RCondThresh',1e-5,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'},'SynapseMemoryModel.Spectrum','RCondThresh'));
% end
% p.parse(varargin{:});
% r=p.Results;



W=obj.GetWf;
q=obj.GetEnc;

[u,qa]=eig(-W);
qa=diag(qa);
[qa,ix]=sort(qa);
u=u(:,ix);

%this method doesn't work when eigenvalues are nearly degenerate
if min(diff(qa)) > obj.DegThresh
    [v,qb]=eig(-W');
    qb=diag(qb);
    [~,ix]=sort(qb);
    v=v(:,ix);
    v=diag(1./diag(v.'*u)) * v.';

    Z=u * diag(1./[1;qa(2:end)]) * v;
    pr=v(1,:);
    pr=pr/sum(pr);
    ca = 2*obj.fp*(1-obj.fp) * (v*obj.w) .* (pr*q*Z*u)';
    return;
end


%this method doesn't work when eigenvectors are nearly parallel
% if rcond(u) > RCondThresh
    Zinv=obj.GetZinv;
    %this method doesn't work when eigenvectors are nearly parallel or W
    %non-ergodic
    if rcond(Zinv) > obj.RCondThresh
        ca = 2*obj.fp*(1-obj.fp) * (u\obj.w) * sum((Zinv\q) * (Zinv\u), 1);
        ca=diag(ca);
    else
        Z=u * diag(1./[1;qa(2:end)]) / u;
        pr=[1 zeros(1,length(qa)-1)] / u;
        pr=pr/sum(pr);
        ca = 2*obj.fp*(1-obj.fp) * (u\obj.w) .* (pr*q*Z*u)';
    end
% else
%     qa=NaN(n,1);
%     ca=qa;
% end
end

