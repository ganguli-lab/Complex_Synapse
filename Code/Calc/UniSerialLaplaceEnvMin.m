function [ smin,smax ] = UniSerialLaplaceEnvMin( numstates )
%[smin,smax]=UNISERIALLAPLACEENVMIN(numstates) extent of maximum Laplace
%transform of SNR curve for uniform serial model
%   Min/max value of s for which envelope computed with UNISERIALLAPLACEENV is valid
%   smin,smax = inverse timescale, Laplace transform parameter
%   numstates = # states

xs = 4.50656136546080;
 smin = ( xs.^(1./numstates) - xs.^(-1./numstates) ).^2 / 2;
 smax = ((xs-1)^2)/(2*xs);

end

