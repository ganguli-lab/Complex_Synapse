function [ genfun ] = MultiExpGenFun( s,data )
%genfun=MULTIEXPGENFUN(s,data) estimate of moment generating function
%   genfun = <exp(-st)>
%   s      = column of s values
%   data   = row of ssamples of t values

genfun = sum(exp(-s*data),2);


end

