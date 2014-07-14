function Validate( obj )
%VALIDATESYNAPSEOPTIMSET(OBJ) validates options in SynapseOptimset OBJ,
%throws error if invalid
%   Parametr name = definition (default,{valid types},{valid attributes})
%     MaxIter = maximum # iterations in fitting (1000,{'numeric'},{'scalar','nonnegative'})
%     TolFun = Stop when log-likelihood reaches this (NaN,{'numeric'},{'scalar'})
%     TolX = Stop when change in model this (1e-4,{'numeric'},{'scalar','nonnegative'})
%     TolFunChange = Stop when change in log-likelihood less than this (1,{'numeric'},{'scalar','nonnegative'})
%     fp = fraction of potentiating events (0.5,{'numeric'},{'row','>=',0,'<=',1})
%     Display = display message at each step? ('off',{'off','iter','final','notify','addstate'})
%     OutputFcn = output function handle ([],{'function_handle},{})
%       syntax: stop=OutputFcn(fitmodelobj,optimValues,state)
%       search for "Output Functions" in the help browser for more information
%     PlotFcn = plot function handle ([],{'function_handle},{})
%       syntax: stop=PlotFcn(fitmodelobj,optimValues,state)
%       search for "Output Functions" in the help browser for more information
%     Algorithm = algorithm used on each sequence ('BW',{'BW','Viterbi'})
%       BW=Baum-Welch, Viterbi=Viterbi path
%     Weighter = algorithm used to combine sequences ('RJ',{'RJ','Uni','Mackay'})
%       RJ=Rabiner-Juang, Uni=uniform weights, Mackay (extra param:
%           Holdback = how many sequences to hold back to evaluate likelihoods
%           positve: which sequences, negative intewger: #sequences,
%           negative fraction: fraction of sequences
%     Penaliser = model penalty (-log prior) ('No',{'No','OffDiagL1','Lhalf'})
%       No=flat prior, OffDiagL1=off-diagonal L^1, Lhalf=L^1/2
%     Penalty = coefficient of model penalty (1,{'numeric'},{'scalar','nonnegative'})
%     ExtraParams = extra param-value pairs passed to algorithms ({},{'cell'},{})
%     ModelDiff = measure of model distance ('KL',{'KL','Ln`'})
%       KL=Kullback-Leibler divergence, Ln=mean L^n norm of rows
%     GroundTruth = true model used to generate data ([],{'SynapseIdModel'},{})
%   For fittng number of states:
%     MaxStates = maximum # states per value of w (6,{'numeric'},{'scalar','nonnegative'})
%     MinLogLikeInc = minimum log likelihood increase to carry on adding states (0,{'numeric'},{'scalar','nonnegative'})
%     NumReps = number of attempts for each w (10,{'numeric'},{'scalar','nonnegative','integer'})

validateattributes(obj.MaxIter,{'numeric'},{'scalar','nonnegative'},'SynapseOptimset','MaxIter');
validateattributes(obj.TolFun,{'numeric'},{'scalar'},'SynapseOptimset','TolFun');
validateattributes(obj.TolX,{'numeric'},{'scalar','nonnegative'},'SynapseOptimset','TolX');
validateattributes(obj.TolFunChange,{'numeric'},{'scalar','nonnegative'},'SynapseOptimset','TolFunChange');
validateattributes(obj.Penalty,{'numeric'},{'scalar','nonnegative'},'SynapseOptimset','Penalty');
validateattributes(obj.fp,{'numeric'},{'row','>=',0,'<=',1},'SynapseOptimset','fp');
validateattributes(obj.ExtraParams,{'cell'},{},'SynapseOptimset','ExtraParams');

if ~isempty(obj.OutputFcn)
    validateattributes(obj.OutputFcn,{'function_handle'},{},'SynapseOptimset','OutputFcn');
end
if ~isempty(obj.PlotFcn)
    validateattributes(obj.PlotFcn,{'function_handle'},{},'SynapseOptimset','PlotFcn');
end
if ~isempty(obj.GroundTruth)
    validateattributes(obj.GroundTruth,{'SynapseIdModel'},{},'SynapseOptimset','GroundTruth');
end
validatestring(obj.Display,{'off','iter','final','notify','addstate'},'SynapseOptimset','Display');
validatestring(obj.Algorithm,{'BW','Viterbi'},'SynapseOptimset','Algorithm');
validatestring(obj.Weighter,{'RJ','Uni','Mackay'},'SynapseOptimset','Weighter');
validatestring(obj.Penaliser,{'No','OffDiagL1','Lhalf'},'SynapseOptimset','Penaliser');
validatestring(obj.ModelDiff,{'KL','Ln`'},'SynapseOptimset','ModelDiff');

validateattributes(obj.MaxStates,{'numeric'},{'scalar','nonnegative','integer'},'SynapseOptimset','MaxStates');
validateattributes(obj.MinLogLikeInc,{'numeric'},{'scalar','nonnegative'},'SynapseOptimset','MinLogLikeInc');
validateattributes(obj.NumReps,{'numeric'},{'scalar','nonnegative','integer'},'SynapseOptimset','NumReps');
validateattributes(obj.NumSample,{'numeric'},{'scalar','nonnegative','integer'},'SynapseOptimset','NumSample');
validateattributes(obj.PriorCcoeff,{'numeric'},{'scalar','nonnegative'},'SynapseOptimset','PriorCcoeff');
validateattributes(obj.PriorQcoeff,{'numeric'},{'scalar','nonnegative'},'SynapseOptimset','PriorQcoeff');


end

