classdef SynapseOptimset < hgsetget
    %options=SYNAPSEOPTIMSET(oldoptions,param,val,...) optimisation parameters for
    %synapse fitting
    %   oldoptions can be omitted, defaults will be used instead
    %   Parametr name = definition (default,{valid types},{valid attributes})
    %                           or (default,{valid vaalues}) for strings
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
    
    properties
        %     MaxIter = maximum # iterations in fitting (1000,{'numeric'},{'scalar','nonnegative'})
        MaxIter=1000;
        %     TolFun = Stop when log-likelihood reaches this (NaN,{'numeric'},{'scalar'})
        TolFun = NaN;
        %     TolX = Stop when change in model this (1e-4,{'numeric'},{'scalar','nonnegative'})
        TolX = 1e-4;
        %     TolFunChange = Stop when change in log-likelihood less than this (1,{'numeric'},{'scalar','nonnegative'})
        TolFunChange = 1;
        %     Display = display message at each step? ('off',{'off','iter','final','notify','addstate'})
        Display = 'off';
        %     OutputFcn = output function handle ([],{'function_handle},{})
        %       syntax: stop=OutputFcn(fitmodelobj,optimValues,state)
        %       search for "Output Functions" in the help browser for more information
        OutputFcn =[];
        %     PlotFcn = plot function handle ([],{'function_handle},{})
        %       syntax: stop=PlotFcn(fitmodelobj,optimValues,state)
        %       search for "Output Functions" in the help browser for more information
        PlotFcn = [];
        %     Algorithm = algorithm used on each sequence ('BW',{'BW','Viterbi'})
        %       BW=Baum-Welch, Viterbi=Viterbi path
        Algorithm = 'BW';
        %     Weighter = algorithm used to combine sequences ('RJ',{'RJ','Uni','Mackay'})
        %       RJ=Rabiner-Juang, Uni=uniform weights, Mackay (extra param:
        %           Holdback = how many sequences to hold back to evaluate likelihoods
        %           positve: which sequences, negative intewger: #sequences,
        %           negative fraction: fraction of sequences
        Weighter = 'RJ';
        %     Penaliser = model penalty (-log prior) ('No',{'No','OffDiagL1','Lhalf'})
        %       No=flat prior, OffDiagL1=off-diagonal L^1, Lhalf=L^1/2
        Penaliser = 'No';
        %     Penalty = coefficient of model penalty (1,{'numeric'},{'scalar','nonnegative'})
        Penalty = 1;
        %     ExtraParams = extra param-value pairs passed to algorithms ({},{'cell'},{})
        ExtraParams = {};
        %     ModelDiff = measure of model distance ('KL',{'KL','Ln`'})
        %       KL=Kullback-Leibler divergence, Ln=mean L^n norm of rows
        ModelDiff = 'Ln';
        %     GroundTruth = true model used to generate data ([],{'SynapseIdModel'},{})
        GroundTruth = [];
        %   For fittng number of states:
        %     MaxStates = maximum # states per value of w (6,{'numeric'},{'scalar','nonnegative','integer'})
        MaxStates = 6;
        %     MinLogLikeInc = minimum log likelihood increase to carry on adding states (0,{'numeric'},{'scalar','nonnegative'})
        %                   or: minimum log likelihood increase on held back data to avoid overfitting
        MinLogLikeInc = 0;
        %     HoldbackForget = forgetting factor for previous likelihood increase on held back data (0,{'numeric'},{'scalar','nonnegative','<',1})
        HoldbackForget = 0;
        %     NumReps = number of attempts for each w (10,{'numeric'},{'scalar','nonnegative','integer'})
        NumReps = 10;
        %     NumSample = number of samples from posterior (10,{'numeric'},{'scalar','nonnegative','integer'})
        NumSample=1000;
        %     PriorCcoeff = coefficient of L^1/2 penalty on exponential
        %     coeffs (6,{'numeric'},{'scalar','nonnegative'})
        PriorCcoeff=50;
        %     PriorQcoeff = coefficient of L^1/2 penalty on exponential
        %     decays (6,{'numeric'},{'scalar','nonnegative'})
        PriorQcoeff=1;
    end
    
    methods
        Validate(obj)
    end
    
    methods%constructor
        function obj=SynapseOptimset(varargin)
            if nargin ~=0%false -> default constructor does nothing
                    %
                    %default parameters:
                    %if we're copying another obj
                    [tempobj,varargin]=extractArgOfType(varargin,'SynapseOptimset');
                    if ~isempty(tempobj)
                        obj=tempobj;
                    end
                    %
                    %set parameter values:
                    %
                    %copy from struct:
                    [Unmatched,varargin]=extractArgOfType(varargin,'struct');
                    if ~isempty(Unmatched)
                        obj=CopyFields(Unmatched,obj);
                    end
                    %
                    %use param/value pairs
                    [obj,varargin]=assignToObject(obj,varargin);
                    %
                    %Extract data:
                    %
            end%if nargin ~=0
        end%function SynapseIdModel
    end%methods constructor

    
end

