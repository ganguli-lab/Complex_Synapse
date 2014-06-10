function [ fitmodel,fval,exitflag,output ] = FitSynapseInit( simobj,guessmodel,options )
%[fitmodel,fval,exitflag,output]=FITSYNAPSEINIT(simobj,guessmodel,options) 
%Fitting Initial of SynapseIdModel (fitmodel) to SynapsePlastSeq simobj,
%without updating M.
%   fval     = log likelihood of fitmodel
%   exitflag = describes the exit condition (1=success, 0=max iterationsm <0=failure) 
%   output   = struct that contains information about the optimization 
%   guessmodel = initial guess SynapseIdModel
%   options    = struct of options

if exist('options','var')
    options=SynapseOptimset(options);
else
    options=SynapseOptimset;
end

%information about the current state of the optimiser
optimValues=struct('iteration',0,'procedure',[options.Algorithm ',' options.Weighter ',' options.Penaliser],'funcCount',0,...
    'NumStates',guessmodel.NumStates,...
    'fval',[],'prev',[],'truth',[],'stateProb',{{}});
state='init';

fitmodel=guessmodel.Sort(options.fp);
optimValues.fval=HMMloglike(fitmodel,simobj);

updater=str2func([options.Penaliser 'penalty']);
extraArgs=[{'Algorithm',['Init' options.Algorithm],'Weighter',options.Weighter},options.ExtraParams];

exitflag=0;
msg=['Exceeded max iterations: ' int2str(options.MaxIter)];

%information about the previous state and current state of the optimiser relative to the previous state
optimValues.prev=struct('model',[],'fval',[],'dfval',[],'KLInitial',[],'KLM',[]);
%information about the current state of the optimiser relative to ground truth
if ~isempty(options.GroundTruth)
    optimValues.truth=struct('model',options.GroundTruth,'fval',[],'dfval',[],'KLInitial',[],'KLM',[]);
    optimValues.truth.fval=HMMloglike(options.GroundTruth,simobj);
    optimValues.truth.dfval=optimValues.truth.fval-optimValues.fval;
end

CallOutFcns;

for i=1:options.MaxIter
    %
    state='interrupt';
    stop=CallOutFcns;
    if stop
        exitflag=-5;
        msg='Stopped by OuputFcn or PlotFcn';
        break;
    end
    %
%store previous model
    %
    optimValues.iteration=i;
    optimValues.prev.model=fitmodel;
    optimValues.prev.fval=optimValues.fval;
    %
%perform the update
    %
    try
    [fitmodel,optimValues.fval,optimValues.stateProb]=updater(optimValues.prev.model,simobj,...
        extraArgs{:});
    fitmodel=fitmodel.Sort(options.fp);
    catch ME
        exitflag=-1;
        msg=['Error: ' ME.message];
        break;
    end
    %
%calculate size of changes in model
    %
    [optimValues.prev.KLInitial,optimValues.prev.KLM]=optimValues.prev.model.KLdivs(fitmodel);
    optimValues.prev.dfval=optimValues.fval-optimValues.prev.fval;
    optimValues.funcCount=i;
    %
%calculate difference of current model from groud truth
    %
    if ~isempty(optimValues.truth)
        %
        optimValues.truth.dfval=optimValues.truth.fval-optimValues.fval;
        [optimValues.truth.KLInitial,optimValues.truth.KLM]=optimValues.truth.model.KLdivs(fitmodel);
        %
%decide whether or not to continue
        %
        if optimValues.truth.dfval < options.TolFun
                if optimValues.truth.KLInitial < options.TolX
                    exitflag=1;
                    msg=['Success. trueloglike - loglike < ' num2str(options.TolFun)...
                        ' and KLdiv from true model to fit model < ' num2str(options.TolX)];
                    break;
                else
                    exitflag=-4;
                    msg=['Not enough data? trueloglike - loglike < ' num2str(options.TolFun)...
                        ' but KLdiv from true Initial to fit Initial > ' num2str(options.TolX)];
                    break;
                end
        end
    end
    %
%decide whether or not to continue
    %
    if optimValues.fval > options.TolFun
        exitflag=1;
        msg=['Success. loglike > ' num2str(options.TolFun)];
        break;
    elseif optimValues.prev.KLInitial < options.TolX && abs(optimValues.prev.dfval) < options.TolFunChange
        if isnan(options.TolFun)
            exitflag=1;
        else
            exitflag=-2;
        end
        msg=['Reached local maximum. Change in loglike < ' num2str(options.TolFunChange)...
            '. KL-div due to change in model < ' num2str(options.TolX)];
        break;
    end
    %
    state='iter';
    stop=CallOutFcns;
    if stop
        exitflag=-5;
        msg='Stopped by OuputFcn or PlotFcn';
        break;
    end
end

fval=optimValues.fval;
output=struct('message',msg,'algortihm',options.Algorithm,'weighter',options.Weighter,'penaliser',options.Penaliser,...
    'iterations',i,'prev',optimValues.prev,'truth',optimValues.truth);
state='done';
CallOutFcns;


    function stopv=CallOutFcns
        %
        stopv=false;
        if ~isempty(options.OutputFcn)
            stopv = options.OutputFcn(fitmodel, optimValues, state);
        end
        %
        if ~isempty(options.PlotFcn)
            stopv = options.PlotFcn(fitmodel, optimValues, state);
        end
    end


end

