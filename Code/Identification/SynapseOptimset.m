function [ options ] = SynapseOptimset( varargin )
%options=SYNAPSEOPTIMSET(oldoptions,param,val,...) Set optimisation parameters for
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
%     NormPower = exponent of L^n norm (2,{'numeric'},{'scalar','nonnegative'})
%     GroundTruth = true model used to generate data ([],{'SynapseIdModel'},{})
%   For fittng number of states:
%     MaxStates = maximum # states per value of w (6,{'numeric'},{'scalar','nonnegative'})
%     MinLogLikeInc = minimum log likelihood increase to carry on adding states (0,{'numeric'},{'scalar','nonnegative'})
%     NumReps = number of attempts for each w (10,{'numeric'},{'scalar','nonnegative'})

options=struct('MaxIter',1000,'TolFun',NaN,'TolX',1e-4,'TolFunChange',1,...
    'Algorithm','BW','Weighter','RJ','Penaliser','No','ExtraParams',{{}},...
    'Penalty',1,'ModelDiff','KL','NormPower',2,...
    'Display','off','OutputFcn',[],'PlotFcn',[],...
    'MaxStates',6,'MinLogLikeInc',0,'NumReps',10,...
    'fp',0.5,'GroundTruth',[]);

if nargin>=1 && isstruct(varargin{1})
    options=UpdateOldOptions(options,varargin{1});
    varargin(1)=[];
end

if length(varargin)>1
    varargin=extractPVpairs(varargin);
    [options,unused]=UpdateOptions(options,varargin);
    options.ExtraParams=[options.ExtraParams unused];
end

ValidateSynapseOptimset(options);

    function [opt,unused]=UpdateOldOptions(defopt,oldopt)
        opt=defopt;
        pr=fields(oldopt);
        unused=cell(1,2*length(pr));
        nextun=1;
        for ii=1:length(pr)
            if ~isempty(oldopt.(pr{ii}))
                if isfield(opt,pr{ii})
                    opt.(pr{ii})=oldopt.(pr{ii});
                else
                    unused{nextun}=pr{ii};
                    unused{nextun+1}=oldopt.(pr{ii});
                    nextun=nextun+2;
                end%if isfield
            end%if ~isempty
        end%for ii
        unused(nextun:end)=[];
    end

    function [opt,unused]=UpdateOptions(oldopt,newoptcell)
        opt=oldopt;
        unused=cell(1,length(newoptcell));
        nextun=1;
        for ii=1:2:length(newoptcell)-1
            if isfield(opt,newoptcell{ii})
                opt.(newoptcell{ii})=newoptcell{ii+1};
            else
                unused{nextun}=newoptcell{ii};
                unused{nextun+1}=newoptcell{ii+1};
                nextun=nextun+2;
            end%if isfield
        end%for ii
        unused(nextun:end)=[];
    end


end

