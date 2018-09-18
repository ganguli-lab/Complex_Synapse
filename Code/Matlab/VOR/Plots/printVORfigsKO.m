function printVORfigsKO( prefix,paramWT,paramKO,df,T_train,T_pre,varargin )
%PRINTVORFIGSKO(prefix,paramWT,paramKO,df,T_train,T_pre) Print all figures
%for WT/KO comparison pre/nopre comparisons, for specified parameter set 
%
%   prefix  = string prepended to all file names, must begin with model name
%   paramWT = parameter used for WT (if different for pot/dep use vector [pot, dep],
%                                   for pooled use [dep_max, dep_min, pot],
%                                   otherwise use scalar for both pot and dep) 
%   paramKO = parameter used for KO dep (for pooled use [dep_max, dep_min])
%   df      = change in fpot for gain inc (if base fpNorm != 0.5, 
%                                         or if different for train/pre,
%                                         use fpots: {fpNorm,fpInc,fpDec})
%   T_train = duration of training
%   T_pre   = duration of pre-training

paramPot = paramWT(1);
param_WT_p = paramWT; %need to keep this for pooled)
paramWT = paramWT(end);
pooled = false;

if strncmpi(prefix,'cascade',3)
    builder_h = @CascadeBuilder;
    n = 10;
elseif strncmpi(prefix,'nonuni',3)
    builder_h = @NonuniBuilder;
    n = 10;
elseif strncmpi(prefix,'serial',3)
    builder_h = @SerialBuilder;
    n = 10;
elseif strncmpi(prefix,'multistate',3)
    builder_h = @MultistateBuilder;
    n = 10;
elseif strncmpi(prefix,'binary',3)
    builder_h = @SerialBuilder;
    n = 2;
elseif strncmpi(prefix,'pooled',3)
    builder_h = @PooledBuilder;
    pooled = true;
    paramPot = param_WT_p(1);
    paramWT = param_WT_p(2:3);
    n = 7;
end

if isscalar(df)
    df = {0.5, 0.5+df, 0.5-df};
end

vexpt=VORbuilderKO(builder_h, n, paramPot, paramWT, paramKO, df{:}, T_train,T_pre, pooled);
vexpt.PrintFigs(prefix);


end

