function [ fitmodel,graphinfo ] = FitSynapseSize( fitsim,testsim,numWvals,varargin )
%FITSYNAPSESIZE Summary of this function goes here
%   Detailed explanation goes here

MaxStates=6;%maximum number of states with each value of w
MinLogLikeInc=0;%carry on adding states if increase of log likelihood is greater than this
NumReps=10;%number of attempts for each w
varargin=assignApplicable(varargin);

NumPlastTypes=max(fitsim.NumPlast,testsim.NumPlast);

graphinfo=struct('numstates',1:(numWvals*MaxStates),'loglike',NaN(1,numWvals*MaxStates));

w=(1:numWvals)';

[fitmodel,loglike]=TryNewW(w);
graphinfo.loglike(length(w))=loglike;
%
disp(length(w));
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
                graphinfo.loglike(length(w))=loglike;
                %
                disp(length(w));
                contadding=true;
            end%if MinLogLikeInc
        end%if MaxStates
    end%for j
end%while

graphinfo.numstates(isnan(graphinfo.loglike))=[];
graphinfo.loglike(isnan(graphinfo.loglike))=[];

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
            [guessmodel,guessloglike]=FitSynapse(fitsim,guessmodel,struct(varargin{:}));
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

