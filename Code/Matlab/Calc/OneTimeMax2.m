function [ s ] = OneTimeMax2( t,t1,n )
%s=ONETIMEMAX2(t,t1,n) curve that maximises SNR at t1, with 2 constraints
%s(t) maximises s(t1).
%   t=full range of times
%   t1=time of max
%   n=#states

if t1<n-1
    s=exp(-t/n-1);
else
    s=(n-1)/t1*exp(-t/t1);
end

end

