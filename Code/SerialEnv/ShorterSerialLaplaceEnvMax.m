function [ smin,smax ] = ShorterSerialLaplaceEnvMax( numstates )
%smax=SHORTERSERIALLAPLACEENVMAX(numstates) extent of maximum Laplace
%transform of SNR curve for uniform serial model
%   Min/max value of s for which envelope computed with STICKYSERIALLAPLACEENV is valid
%   smin,smax = inverse timescale, Laplace transform parameter
%   numstates = # states


xs = 4.50656136546080;

alphamin=fzero(@qmax1,xs^(2/numstates));

smin=(alphamin-1)^2/(2*alphamin);

alphamax=fzero(@qmax0,xs^(2/(numstates-2)));

smax=(alphamax-1)^2/(2*alphamax);

    function f=qmax1(alpha)
        f=(-1).*((-1)+alpha.^((1/2).*numstates)).^2.*(1+ ...
          alpha.^numstates)+((-1)+alpha).*alpha.^((-2)+(1/2).* ...
          numstates).*(1+((-1)+alpha).*alpha).*((-1).*alpha+ ...
          alpha.^numstates).*numstates;
    end

    function f=qmax0(alpha)
        f=(-1)+alpha+(-1).*alpha.^2+2.*alpha.^(1+(1/2).*numstates)+( ...
  -2).*alpha.^((-2)+numstates)+2.*alpha.^((-1)+numstates)+(-2) ...
  .*alpha.^numstates+2.*alpha.^((-3)+(3/2).*numstates)+(-1).* ...
  alpha.^((-4)+2.*numstates)+alpha.^((-3)+2.*numstates)+(-1).* ...
  alpha.^((-2)+2.*numstates)+alpha.^((-1)+(1/2).*numstates).* ...
  numstates+(-1).*alpha.^((1/2).*numstates).*numstates+(-1).* ...
  alpha.^((-2)+(3/2).*numstates).*numstates+alpha.^((-1)+(3/2) ...
  .*numstates).*numstates;
    end


end

