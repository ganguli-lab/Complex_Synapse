function [ comps ] = ScanCNNum( ranges, ns, reps, useCnotN, varargin )
%comps=SCANPOOLEDNUM(ranges,ns) parameter scan for pooled resource model
%   comps  = learning rate differences: WT_nopre - WT_pre
%   ranges = range for each parameter scan
%   ns     = list of numbers of states
%   reps   = number of attempts
%   parametrs: pot, dep_min, dep_max, fp_norm, fp_inc, fp_dec

comps = -Inf(2,length(ns));
DispCounter(1,ns(end),'n:');
for n = ns
    DispCounter(n,ns(end),'n:');
    pcomps = ScanCNtop(ranges, n, reps, useCnotN, varargin{:});
%     pcomps = ScanCNtop2(ranges, n, useCnotN, varargin{:});
    comps(1,ns==n) = max(pcomps(~isnan(pcomps)));
    comps(2,ns==n) = min(pcomps(~(isnan(pcomps) | isinf(pcomps))));
end
DispCounter(ns(end)+1,ns(end),'n:');
end

