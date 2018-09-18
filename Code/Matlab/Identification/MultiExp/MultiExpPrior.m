function [ neglogprior ] = MultiExpPrior( params,penc,penq )
%neglogprior=MULTIEXPPRIOR(arams,penc,penq) prior distribution on
%parameters of multi-exponential distribution
%   params = column [log(q); c(1:end-1)]
%       c(end)=1-sum(c(1:end-1))
%       like = prod_data sum_{c,q} c*q*exp(-q*data)
%   prior = exp(-penc*||c||_1/2) exp(-penq*||q||_2)
%   neglogprior = -log(prior)

% persistent p
% if isempty(p)
%     p=inputParser;
%     p.FunctionName='MultiExpPrior';
%     p.StructExpand=true;
%     p.KeepUnmatched=false;
%     p.addOptional('penc',1,...
%         @(x)validateattributes(x,{'numeric'},{'scalar'},'MultiExpPrior','penc',2));
%     p.addOptional('penq',1,...
%         @(x)validateattributes(x,{'numeric'},{'scalar'},'MultiExpPrior','optimOptions',3));
% end
% p.parse(varargin{:});
% r=p.Results;

n=ceil(length(params)/2);
q2=exp(2*params(1:n));
% tau=exp(-params(1:n));
c=params(n+1:end);
c=[c;1-sum(c)];

neglogprior = penc*sum(sqrt(abs(c))) + sum(penq*q2/2 - params(1:n));

if min(c) < 0
    neglogprior = inf;
end
% for i = 1:n
%     t_test=linspace(tau(i)/10,tau(i)*10,100);
%     prob=MultiExpLike(params,t_test);
%     if isinf(prob)
%         neglogprior=inf;
%     end
% end

end

