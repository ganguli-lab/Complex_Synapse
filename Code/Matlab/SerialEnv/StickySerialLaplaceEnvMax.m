function [ smax ] = StickySerialLaplaceEnvMax( numstates )
%smax=STICKYSERIALLAPLACEENVMAX(numstates) extent of maximum Laplace
%transform of SNR curve for uniform serial model
%   Min/max value of s for which envelope computed with STICKYSERIALLAPLACEENV is valid
%   smax      = inverse timescale, Laplace transform parameter
%   numstates = # states


xs = 4.50656136546080;

alphamax=fzero(@dummy,xs^(2/numstates));

smax=(alphamax-1)^2/(2*alphamax);

    function f=dummy(alpha)
        f=1+(-2).*alpha.^((1/2).*numstates)+2.*alpha.^numstates+(-2).* ...
          alpha.^((3/2).*numstates)+alpha.^(2.*numstates)+alpha.^(1+( ...
          1/2).*numstates).*numstates+(-1).*alpha.^((1/2).*numstates) ...
          .*numstates+(-1).*alpha.^((3/2).*numstates).*numstates+ ...
          alpha.^((-1)+(3/2).*numstates).*numstates;
    end


end

