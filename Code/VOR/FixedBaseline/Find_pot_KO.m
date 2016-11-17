function [ pot_KO ] = Find_pot_KO( builder_h, n, pot_WT, dep_WT, dep_KO, fpNorm, varargin )
%pot_KO=FIND_POT_KO(builder_h,n,pot_WT,dep_WT,dep_KO,fpNorm,reps,minval,maxval)
%Find the value of pot_KO that results in WT and KO having the same
%baseline mean synaptic weight
%   builder_h = function handle that builds Wp, Wm and w
%   n      = number of states
%   pot_?? = parameter for potentiation
%   dep_?? = parameter for depression
%   ???_WT = wild-type parameter
%   ???_KO = wild-type parameter
%   fpNorm = fraction of potentiation events at baseline
%   reps   = number of attempts
%   m??val = min/max value for pot_KO

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='Find_pot_KO';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addOptional('reps',20,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','integer'},'Find_pot_KO','reps'));
    p.addOptional('minval',0,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'},'Find_pot_KO','minval'));
    p.addOptional('maxval',1,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'},'Find_pot_KO','maxval'));
end
p.parse(varargin{:});
r=p.Results;

bw = BaselineWt(builder_h, n, pot_WT, dep_WT, fpNorm);
lossfun = @(x) (BaselineWt(builder_h, n, x, dep_KO, fpNorm) - bw)^2;

opts = optimoptions('fmincon','Display','off');

x0 = r.minval + (r.maxval - r.minval) * rand;
lb = r.minval; 
ub = r.maxval; 

minf = Inf;
pot_KO = NaN;

for i = 1:r.reps
     [x, fval, ef] = fmincon(lossfun, x0, [],[],[],[], lb, ub, [], opts);
%     [x, fval] = fmincon(lossfun, x0, [],[],[],[], lb, ub, [], opts);
    if fval < minf && fval < 1e-3 && ef > 0 
        pot_KO = x;
        minf = fval;
    end
end

if ~isnan(pot_KO)
    y = BaselineWt(builder_h, n, pot_KO, dep_KO, fpNorm);
    if (y - bw)^2 >= 1e-3 || pot_KO < lb || pot_KO > ub
        error('invalid result:');
    end
end


end

