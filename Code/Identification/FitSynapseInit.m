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

[optimValues,updaterfn,extraArgs]=FitSynapseUtility_setup(options,guessmodel,simobj,'Init');

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
    [optimValues,exitflag,msg,ME]=FitSynapseUtility_update(optimValues,options,simobj,updaterfn,extraArgs);
    %
    CallOutFcns;
    if exitflag~=0
        break;
    end
    %
    optimValues=FitSynapseUtility_calcchanges(optimValues,options);
    %
    CallOutFcns;
    if exitflag~=0
        break;
    end
    %
    [optimValues,exitflag,msg]=FitSynapseUtility_checktermination(optimValues,options,'Init');
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

[fitmodel,fval,output]=FitSynapseUtility_finish(optimValues,options,msg,ME);

state='done';
CallOutFcns;


    function stopv=CallOutFcns
        %
        stopv=false;
        if ~isempty(options.OutputFcn)
            stopv = options.OutputFcn(optimValues.model, optimValues, state);
        end
        %
        if ~isempty(options.PlotFcn)
            stopv = options.PlotFcn(optimValues.model, optimValues, state);
        end
        if stopv
            exitflag=-5;
            msg='Stopped by OuputFcn or PlotFcn';
        end
        FitSynapseUtility_display(optimValues,options,state,exitflag,msg);
    end


end

