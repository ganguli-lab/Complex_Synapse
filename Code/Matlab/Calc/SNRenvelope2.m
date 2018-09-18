function [ s,bnds ] = SNRenvelope2( t,n )
%[s,bnds]=SNRENVELOPE(t,n) maximaum possible SNR at each time, with two
%constraints
%   t = time (row vector)
%   n = number of states
%   s = SNR enveope
%   bnds = boundaries of regions where different constraints are active, column vector. 


bnds=n-1;

s1=exp(-t/(n-1));
s2=(n-1)./(exp(1)*t);

t1= t<=bnds(1);
t2= t>bnds(1);

s=s1.*t1 + s2.*t2;

end

