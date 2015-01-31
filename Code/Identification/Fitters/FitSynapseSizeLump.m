function [ fitmodel ] = FitSynapseSizeLump( simobjs,varargin )
%[fitmodel,like_n]=FitSynapseSizeLump(simobjs,options)
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
    p.KeepUnmatched=true;
    p.addOptional('options',SynapseOptimset,@(x)validateattributes(x,{'SynapseOptimset'},{},'FitSynapseSize','options',2));
end
p.parse(varargin{:});
options=p.Results.options;

NumPlastTypes=simobjs.NumPlast;
numWvals=simobjs.NumWvals;

w=zeros(numWvals*options.MaxStates,1);
w(1:options.MaxStates:end)=1;
w=cumsum(w);

fitmodel=SynapseIdModel.Rand(w,'NumPlastTypes',NumPlastTypes,p.Unmatched);
testloglike=HMMloglike(fitmodel,simobjs)+SynapsePrior(fitmodel,options);
for i=1:options.NumReps
    DispReps(i);
    guessmodel=SynapseIdModel.Rand(w,'NumPlastTypes',NumPlastTypes,p.Unmatched);
    [guessmodel,guessloglike]=FitSynapse(simobjs,guessmodel,options);
    if guessloglike > testloglike
        fitmodel=guessmodel;
        testloglike=guessloglike;
    end%if guessloglike
end%for i
DispReps(i+1);

partitions=fitmodel.FindLumps;
if fitmodel.TestLump(partitions);
    fitmodel=fitmodel.Lumpify(partitions);
end




    function DispReps(i)
        if ~strcmpi(options.Display,'off')
            if i>1 && strcmpi(options.Display,'addstate')
                fprintf(repmat('\b',[1 numel(sprintf('Rep:%d/%d ',i-1,options.NumReps))]));
            end
            if i<=options.NumReps
                fprintf('Rep:%d/%d\n',i,options.NumReps)
            end
        end
    end



end

