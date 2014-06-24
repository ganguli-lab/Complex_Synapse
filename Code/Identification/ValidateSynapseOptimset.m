function ValidateSynapseOptimset( options )
%VALIDATESYNAPSEOPTIMSET(options) validates options from SynapseOptimset,
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
%     NumReps = number of attempts for each w (10,{'numeric'},{'scalar','nonnegative'})

validateattributes(options.MaxIter,{'numeric'},{'scalar','nonnegative'},'SynapseOptimset','MaxIter');
validateattributes(options.TolFun,{'numeric'},{'scalar'},'SynapseOptimset','TolFun');
validateattributes(options.TolX,{'numeric'},{'scalar','nonnegative'},'SynapseOptimset','TolX');
validateattributes(options.TolFunChange,{'numeric'},{'scalar','nonnegative'},'SynapseOptimset','TolFunChange');
validateattributes(options.Penalty,{'numeric'},{'scalar','nonnegative'},'SynapseOptimset','Penalty');
validateattributes(options.MaxStates,{'numeric'},{'scalar','nonnegative'},'SynapseOptimset','MaxStates');
validateattributes(options.MinLogLikeInc,{'numeric'},{'scalar','nonnegative'},'SynapseOptimset','MinLogLikeInc');
validateattributes(options.NumReps,{'numeric'},{'scalar','nonnegative'},'SynapseOptimset','NumReps');
validateattributes(options.fp,{'numeric'},{'row','>=',0,'<=',1},'SynapseOptimset','fp');
validateattributes(options.ExtraParams,{'cell'},{},'SynapseOptimset','ExtraParams');
if ~isempty(options.OutputFcn)
    validateattributes(options.OutputFcn,{'function_handle'},{},'SynapseOptimset','OutputFcn');
end
if ~isempty(options.PlotFcn)
    validateattributes(options.PlotFcn,{'function_handle'},{},'SynapseOptimset','PlotFcn');
end
if ~isempty(options.GroundTruth)
    validateattributes(options.GroundTruth,{'SynapseIdModel'},{},'SynapseOptimset','GroundTruth');
end
validatestring(options.Display,{'off','iter','final','notify','addstate'},'SynapseOptimset','Display');
validatestring(options.Algorithm,{'BW','Viterbi'},'SynapseOptimset','Algorithm');
validatestring(options.Weighter,{'RJ','Uni','Mackay'},'SynapseOptimset','Weighter');
validatestring(options.Penaliser,{'No','OffDiagL1','Lhalf'},'SynapseOptimset','Penaliser');
validatestring(options.ModelDiff,{'KL','Ln`'},'SynapseOptimset','ModelDiff');


end

