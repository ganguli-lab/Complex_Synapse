function [ fitmodel,like_n ] = FitSynapseSize( fitsim,testsim,varargin )
%[fitmodel,like_n]=FitSynapseSize(fitsim,testsim,varargin)
%Determinig number of states needed for SynapseIdModel (fitmodel) to fit
%SypssePlastSeq (fitsim,testsim)
%   like_n = struct(numstates,loglike)
%   fitsim  = used to fit fitmodel.M
%   testsim = used to evaluate models, after refitting fitmodel.Initial

MaxStates=6;%maximum number of states with each value of w
MinLogLikeInc=chi2inv(0.95,1)/2;%carry on adding states if increase of log likelihood is greater than this
NumReps=10;%number of attempts for each w
Display=true;
varargin=assignApplicable(varargin);

NumPlastTypes=max(fitsim.NumPlast,testsim.NumPlast);
numWvals=max(fitsim.NumWvals,testsim.NumWvals);

like_n=struct('numstates',1:(numWvals*MaxStates),'loglike',NaN(1,numWvals*MaxStates));

w=(1:numWvals)';

[fitmodel,loglike]=TryNewW(w);
like_n.loglike(length(w))=loglike;
%
if Display
    disp(length(w));
end
contadding=true;

while contadding
    contadding=false;
    for j=1:numWvals
        if sum(w==j) < MaxStates
            neww=AddState(w,j);
%             newmodel=TryNewW(neww);
            [newmodel,newloglike]=TryNewW(neww);
%             newloglike=HMMloglike(newmodel,testsim);
            if newloglike-loglike > MinLogLikeInc
                w=neww;
                fitmodel=newmodel;
                loglike=newloglike;
                like_n.loglike(length(w))=loglike;
                %
                if Display
                    disp(length(w));
                end
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
        newmodel=FitModel(neww);
        [newmodel,newloglike]=ReFitInitial(newmodel);
    end%function TryNewW

    function newmodel=FitModel(neww)
        newmodel=SynapseIdModel.Rand(neww);
        testloglike=HMMloglike(newmodel,fitsim);
        for i=1:NumReps
            guessmodel=SynapseIdModel.Rand(neww,'NumPlastTypes',NumPlastTypes);
            [guessmodel,guessloglike]=FitSynapseM(fitsim,guessmodel,struct(varargin{:}));
            if guessloglike > testloglike
                newmodel=guessmodel;
                testloglike=guessloglike;
            end%if guessloglike
        end%for i
    end%function TryNewW

    function [newmodel,newloglike]=ReFitInitial(trialmodel)
        newloglike=-Inf;
        newmodel=trialmodel;
        for i=1:NumReps
            p=rand(1,trialmodel.NumStates);
            guessmodel=trialmodel.setInitial(p/sum(p));
            [guessmodel,guessloglike]=FitSynapseInit(testsim,guessmodel,struct(varargin{:}));
            if guessloglike > newloglike
                newmodel=guessmodel;
                newloglike=guessloglike;
            end%if guessloglike
        end%for i
    end%function TryNewW



end

