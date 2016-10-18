function [ comps ] = ScanPooledNum( ranges, ns )
%comps=SCANPOOLEDLEFT(ranges,n) parameter scan for pooled resource model
%   comps  = learning rate differences: WT_nopre - WT_pre
%   ranges = range for each parameter scan
%   ns     = list of numbers of states
%   parametrs: pot, dep_min, dep_max, fp_norm, fp_inc, fp_dec

comps = -Inf(size(ns));
DispCounter(1,ns(end),'n:');
for n = ns
    DispCounter(n,ns(end),'n:');
    pcomps = ScanPooledLeft(ranges, n);
    comps(ns==n) = max(pcomps(:));
end
DispCounter(ns(end)+1,ns(end),'n:');
end

