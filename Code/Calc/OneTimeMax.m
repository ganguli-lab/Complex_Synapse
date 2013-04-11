function [ s ] = OneTimeMax( t,t1,n )
%s=ONETIMEMAX(t,t1,n) curve that maximises SNR at t1, with 3 constraints
%s(t) maximises s(t1).
%   t=full range of times
%   t1=time of max
%   n=#states
gammasq=(128)/pi^4;

if t1<gammasq/2
    s=exp(-t/gammasq);
elseif t1<(n-1)^2/(2*gammasq)
    s=sqrt(gammasq/(2*t1))*exp(-t/(2*t1));
elseif t1<(n-1)^2/(gammasq)
    s=gammasq/(n-1)*exp(-t*gammasq/(n-1)^2);
else
    s=(n-1)/t1*exp(-t/t1);
end



end

