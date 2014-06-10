function [ fitmodel,fval,exitflag,output ] = FitSynapse( simobj,guessmodel,options )
%[fitmodel,fval,exitflag,output]=FITSYNAPSE(simobj,guessmodel,options)
%Fitting a SynapseIdModel (fitmodel) to SynapsePlastSeq simobj
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

[optimValues,fitmodel,updaterfn,extraArgs]=FitSynapseUtility_setup(options,guessmodel,simobj);

exitflag=0;
msg=['Exceeded max iterations: ' int2str(options.MaxIter)];

state='init';
CallOutFcns;

while exitflag==0 &&  optimValues.iteration <= options.MaxIter
    %
    state='interrupt';
    CallOutFcns;
    if exitflag~=0
        break;
    end
    %
    [optimValues,fitmodel,exitflag,msg]=FitSynapseUtility_update(optimValues,options,fitmodel,simobj,updaterfn,extraArgs);
    %
    CallOutFcns;
    if exitflag~=0
        break;
    end
    %
    optimValues=FitSynapseUtility_calcchanges(optimValues,fitmodel);
    %
    CallOutFcns;
    if exitflag~=0
        break;
    end
    %
    [optimValues,exitflag,msg]=FitSynapseUtility_checktermination(optimValues,options);
    %
    CallOutFcns;
    if exitflag~=0
        break;
    end
    %
    state='iter';
    CallOutFcns;
    %
end%while

[fval,output]=FitSynapseUtility_finish(optimValues,options,msg);

state='done';
CallOutFcns;


    function CallOutFcns
        %
        stopv=false;
        if ~isempty(options.OutputFcn)
            stopv = options.OutputFcn(fitmodel, optimValues, state);
        end
        %
        if ~isempty(options.PlotFcn)
            stopv = options.PlotFcn(fitmodel, optimValues, state);
        end
        if stopv
            exitflag=-5;
            msg='Stopped by OuputFcn or PlotFcn';
        end
    end


end

