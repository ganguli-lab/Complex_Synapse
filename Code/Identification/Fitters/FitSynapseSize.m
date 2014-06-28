function [ fitmodel,like_n ] = FitSynapseSize( simobjs,options,varargin )
%[fitmodel,like_n]=FitSynapseSize(simobjs,options)
%Determinig number of states needed for SynapseIdModel (fitmodel) to fit
%SypssePlastSeq (simobjs)
%   like_n = struct(numstates,loglike)
%   fitsim  = used to fit fitmodel.M
%   testsim = used to evaluate models, after refitting fitmodel.Initial
%other arguments passed to SynapseIdModel.Rand

% MaxStates=6;%maximum number of states with each value of w
% MinLogLikeInc=0;%carry on adding states if increase of log likelihood is greater than this
% MinLogLikeInc=chi2inv(0.95,1)/2;%carry on adding states if increase of log likelihood is greater than this
% NumReps=10;%number of attempts for each w
% Display=true;
% varargin=assignApplicable(varargin);
if exist('options','var')
    options=SynapseOptimset(options);
else
    options=SynapseOptimset;
end

NumPlastTypes=simobjs.NumPlast;
numWvals=simobjs.NumWvals;

like_n=struct('numstates',1:(numWvals*options.MaxStates),'loglike',NaN(1,numWvals*options.MaxStates));

w=(1:numWvals)';

[fitmodel,loglike]=TryNewW(w);
like_n.loglike(length(w))=loglike;
%
DispStates;
contadding=true;

while contadding
    contadding=false;
    for j=1:numWvals
        if sum(w==j) < options.MaxStates
            neww=AddState(w,j);
            [newmodel,newloglike]=TryNewW(neww);
            if newloglike-loglike > options.MinLogLikeInc
                w=neww;
                fitmodel=newmodel;
                loglike=newloglike;
                like_n.loglike(length(w))=loglike;
                %
                DispStates;
                contadding=true;
            end%if MinLogLikeInc
        end%if MaxStates
    end%for j
end%while

like_n.numstates(isnan(like_n.loglike))=[];
like_n.loglike(isnan(like_n.loglike))=[];

    function neww=AddState(oldw,wval)
        i=find(oldw==wval,1,'first');
        neww=[oldw(1:i-1); wval; oldw(i:end)];
    end%function addstate

    function [newmodel,newloglike]=TryNewW(neww)
        newloglike=zeros(size(simobjs,1),1);
        for k=1:size(simobjs,1);
            inds=1:length(newloglike);
            inds(k)=[];
            newmodel=FitModel(neww,reshape(simobjs(inds,:),1,[]));
            [newmodel,newloglike(k)]=ReFitInitial(newmodel,simobjs(k,:));
        end%for k
        newloglike=mean(newloglike);
    end%function TryNewW

    function newmodel=FitModel(neww,fitsim)
        newmodel=SynapseIdModel.Rand(neww);
        testloglike=HMMloglike(newmodel,fitsim)+SynapsePrior(newmodel,options);
        for i=1:options.NumReps
            guessmodel=SynapseIdModel.Rand(neww,'NumPlastTypes',NumPlastTypes,varargin{:});
            [guessmodel,guessloglike]=FitSynapseM(fitsim,guessmodel,options);
            if guessloglike > testloglike
                newmodel=guessmodel;
                testloglike=guessloglike;
            end%if guessloglike
        end%for i
    end%function TryNewW

    function [newmodel,newloglike]=ReFitInitial(trialmodel,testsim)
        newloglike=-Inf;
        newmodel=trialmodel;
        for i=1:options.NumReps
            p=rand(1,trialmodel.NumStates);
            guessmodel=trialmodel.setInitial(p/sum(p));
            [guessmodel,guessloglike]=FitSynapseInit(testsim,guessmodel,options);
            if guessloglike > newloglike
                newmodel=guessmodel;
                newloglike=guessloglike;
            end%if guessloglike
        end%for i
    end%function TryNewW

    function DispStates
        if ~strcmpi(options.Display,'off')
            disp(['# states: ' int2str(length(w))]);
        end
    end


end

