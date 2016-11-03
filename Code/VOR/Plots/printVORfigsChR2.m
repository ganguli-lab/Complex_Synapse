function printVORfigsChR2( prefix,param,df,T_train,T_pre,varargin )
%PRINTVORFIGSCHR2(prefix,param,df,T_train,T_pre) Print all figures
%for CF stim/no CF stim comparisons, for specified parameter set 
%
%   prefix  = string prepended to all file names, must begin with model name
%   param   = parameter used for WT (if different for pot/dep use vector [pot, dep],
%                                   for pooled use [dep_max, dep_min, pot],
%                                   otherwise use scalar for both pot and dep) 
%   df      = change in fpot for gain inc (if base fpNorm != 0.5, 
%                                         or if different for train/pre,
%                                         use fpots: {fpNorm,fpInc,fpDec})
%   T_train = duration of training
%   T_pre   = duration of pre-training

paramPot = param(1);
param_p = param; %need to keep this for pooled
paramDep = param(end);
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
    paramPot = param_p(1);
    paramDep = param_p(2:3);
    n = 7;
end

if isscalar(df)
    df = {0.5, 0.5+df, (0.5+df)^2};
end

vexpt=VORbuilderChR2(builder_h, n, paramPot, paramDep, df{:}, T_train, T_pre, pooled);
% vexpt=VORbuilderChR2eq(builder_h, n, paramPot, paramDep, df{:}, T_train, pooled);

vexpt.LegFontSize=2*vexpt.LegFontSize;
vexpt.LabFontSize=2*vexpt.LabFontSize;
vexpt.txFontSize=2*vexpt.txFontSize;
vexpt.EqFontSize=2*vexpt.EqFontSize;
vexpt.ProbFontSize=4*vexpt.ProbFontSize;


fig=figure('PaperPositionMode','auto','Position',[60 60 1000 1000]);
Parent=axes('Parent',fig);
vexpt.PlotLearn('LineWidth',2,'Parent',Parent);
print(fig,[prefix '_ChR2_learn.eps'],'-depsc');
close(fig);

figs=figure('PaperPositionMode','auto','Position',[60 60 1000 1000]);
Parent=axes('Parent',figs);
vexpt.PlotLearnS('LineWidth',2,'Parent',Parent);
print(figs,[prefix '_ChR2_learnS.eps'],'-depsc');
close(figs);
    
    
end

