function [ fitmodel,like_n ] = FitSynapseSize( simobjs,varargin )
%[fitmodel,like_n]=FitSynapseSize(simobjs,options)
%Determinig number of states needed for SynapseIdModel (fitmodel) to fit
%SypssePlastSeq (simobjs)
%   like_n = struct(numstates,loglike)
%   fitsim  = used to fit fitmodel.M
%   testsim = used to evaluate models, after refitting fitmodel.Initial
%other arguments passed to SynapseIdModel.Rand

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='FitSynapseSize';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addOptional('options',SynapseOptimset,@(x)validateattributes(x,{'SynapseOptimset'},{},'FitSynapseSize','options',2));
    p.addOptional('extraArgs',{},@(x)validateattributes(x,{'cell'},{},'FitSynapseSize','extraArgs',3));
%     p.addParameter('Normalise',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
end
p.parse(varargin{:});
options=p.Results.options;

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
            guessmodel=SynapseIdModel.Rand(neww,'NumPlastTypes',NumPlastTypes,p.Results.extraArgs{:});
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

