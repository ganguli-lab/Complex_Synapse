function [ s,bnds ] = SNRenvelope( t,n )
%[s,bnds]=SNRENVELOPE(t,n) maximaum possible SNR at each time
%   t = time (row vector)
%   n = number of states
%   s = SNR enveope
%   bnds = boundaries of regions where different constraints are active, column vector. 

gamma=sqrt(128)/pi^2;

bnds=[gamma^2/2;(n-1)^2/(2*gamma^2);(n-1)^2/(gamma^2);];

s1=exp(-t/gamma^2);
s2=gamma./sqrt(2*exp(1)*t);
s3=gamma^2/(n-1)*exp(-gamma^2*t/(n-1)^2);
s4=(n-1)./(exp(1)*t);

t1= t<=bnds(1);
t2= t>bnds(1) & t<=bnds(2);
t3= t>bnds(2) & t<=bnds(3);
t4= t>bnds(3);

s=s1.*t1 + s2.*t2 + s3.*t3 + s4.*t4;

end

