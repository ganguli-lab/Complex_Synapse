function [ fitmodel,fval,exitflag,output ] = FitSynapseInit( simobj,guessmodel,options )
%[fitmodel,fval,exitflag,output]=FITSYNAPSEINIT(simobj,guessmodel,options) 
%Fitting Initial of SynapseIdModel (fitmodel) to SynapsePlastSeq simobj,
%without updating M.
%   fval     = log likelihood of fitmodel
%   exitflag = describes the exit condition (1=success, 0=max iterationsm <0=failure) 
%   output   = struct that contains information about the optimization 
%   guessmodel = initial guess SynapseIdModel
%   options    = struct of options

defoptions=struct('MaxIter',1000,'TolFun',NaN,'TolX',1e-4,'TolFunChange',1,...
    'Algorithm','BW','Weighter','RJ','ExtraParams',{{}},...
    'Display','off','OutputFcn',[],'PlotFcn',[],...
    'fp',0.5,'GroundTruth',[]);
if exist('options','var')
    [defoptions,unused]=UpdateOptions(defoptions,options);
    defoptions.ExtraParams=[defoptions.ExtraParams unused];
end

optimValues=struct('iteration',0,'procedure',[defoptions.Algorithm ',' defoptions.Weighter],'funcCount',0,...
    'NumStates',guessmodel.NumStates,...
    'fval',[],'prev',[],'truth',[],'stateProb',{{}});
state='init';

fitmodel=guessmodel.Sort(defoptions.fp);
optimValues.fval=HMMloglike(fitmodel,simobj);

weighter=str2func([defoptions.Weighter 'weight']);

exitflag=0;
msg=['Exceeded max iterations: ' int2str(defoptions.MaxIter)];

optimValues.prev=struct('model',[],'fval',[],'dfval',[],'KLInitial',[],'KLM',[]);
if ~isempty(defoptions.GroundTruth)
    optimValues.truth=struct('model',defoptions.GroundTruth,'fval',[],'dfval',[],'KLInitial',[],'KLM',[]);
    optimValues.truth.fval=HMMloglike(defoptions.GroundTruth,simobj);
    optimValues.truth.dfval=optimValues.truth.fval-optimValues.fval;
end

CallOutFcns;

for i=1:defoptions.MaxIter
    %
    state='interrupt';
    stop=CallOutFcns;
    if stop
        exitflag=-5;
        msg='Stopped by OuputFcn or PlotFcn';
        break;
    end
    %
    optimValues.iteration=i;
    optimValues.prev.model=fitmodel;
    optimValues.prev.fval=optimValues.fval;
    %
    [fitmodel,optimValues.fval,optimValues.stateProb]=weighter(optimValues.prev.model,simobj,...
        'Algorithm',['Init' defoptions.Algorithm],defoptions.ExtraParams{:});
    fitmodel=fitmodel.Sort(defoptions.fp);
    %
    [optimValues.prev.KLInitial,optimValues.prev.KLM]=optimValues.prev.model.KLdivs(fitmodel);
    optimValues.prev.dfval=optimValues.fval-optimValues.prev.fval;
    optimValues.funcCount=i;
    %
    if ~isempty(optimValues.truth)
        %
        optimValues.truth.dfval=optimValues.truth.fval-optimValues.fval;
        [optimValues.truth.KLInitial,optimValues.truth.KLM]=optimValues.truth.model.KLdivs(fitmodel);
        %
        if optimValues.truth.dfval < defoptions.TolFun
            if mean(optimValues.truth.KLM/optimValues.NumStates) < defoptions.TolX
                if optimValues.truth.KLInitial < defoptions.TolX
                    exitflag=1;
                    msg=['Success. trueloglike - loglike < ' num2str(defoptions.TolFun)...
                        ' and KLdiv from true model to fit model < ' num2str(defoptions.TolX)];
                    break;
                else
                    exitflag=-4;
                    msg=['Not enough data? trueloglike - loglike < ' num2str(defoptions.TolFun)...
                        ' and KLdiv from true M to fit M > ' num2str(defoptions.TolX)...
                        ' but not Initial'];
                    break;
                end
            else
                exitflag=-3;
                msg=['Not enough data? trueloglike - loglike < ' num2str(defoptions.TolFun)...
                    ' despite KLdiv from true M to fit M > ' num2str(defoptions.TolX)];
                break;
            end
        end
    end
    %
    if optimValues.fval > defoptions.TolFun
        exitflag=1;
        msg=['Success. loglike > ' num2str(defoptions.TolFun)];
        break;
    elseif optimValues.prev.KLInitial < defoptions.TolX && abs(optimValues.prev.dfval) < defoptions.TolFunChange
        if isnan(defoptions.TolFun)
            exitflag=1;
        else
            exitflag=-2;
        end
        msg=['Reached local maximum. Change in loglike < ' num2str(defoptions.TolFunChange)...
            '. KL-div due to change in model < ' num2str(defoptions.TolX)];
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
output=struct('message',msg,'algortihm',defoptions.Algorithm,'weighter',defoptions.Weighter,...
    'iterations',i,'prev',optimValues.prev,'truth',optimValues.truth);
state='done';
CallOutFcns;


    function [opt,unused]=UpdateOptions(oldopt,newopt)
        opt=oldopt;
        pr=fields(newopt);
        unused=cell(1,2*length(pr));
        nextun=1;
        for ii=1:length(pr)
            if ~isempty(newopt.(pr{ii}))
                if isfield(opt,pr{ii})
                    opt.(pr{ii})=newopt.(pr{ii});
                else
                    unused{nextun}=pr(ii);
                    unused{nextun+1}=newopt.(pr{ii});
                    nextun=nextun+2;
                end%if isfield
            end%if ~isempty
        end%for ii
        unused(nextun:end)=[];
    end

    function stopv=CallOutFcns
        %
        stopv=false;
        if ~isempty(defoptions.OutputFcn)
            stopv = defoptions.OutputFcn(fitmodel, optimValues, state);
        end
        %
        if ~isempty(defoptions.PlotFcn)
            stopv = defoptions.PlotFcn(fitmodel, optimValues, state);
        end
    end


end

