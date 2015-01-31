function [ optimValues,exitflag,msg ] = FitSynapseUtility_checktermination( optimValues,options,fitMorInit )
%[optimValues,exitflag,msg]=FITSYNAPSEUTILITY_CHECKTERMINATION(optimValues,options,fitMorInit)
%check if optimiser should terminate
%   optimValues = struct with information about the current state of the optimiser
%   options     = struct of options (see SynapseOptimset)
%   fitMorInit  = 'M' if we're only fitting SynapseIdModel.M, or 'Init' if
%   we're only fitting SynapseIdModel.Initial 
%   exitflag    = integer describing reason for exiting
%   msg         = string describing reason for exiting

fitM=true;
fitInit=true;
if exist('fitMorInit','var')
    switch fitMorInit
     case 'M'
         fitInit=false;
     case 'Init'
         fitM=false;
    end%switch fitMorInit
end%if exist fitMorInit

exitflag=0;
msg=['Exceeded max iterations: ' int2str(options.MaxIter)];
    %
    if ~isempty(optimValues.truth)
        %
%decide whether or not to continue
%using difference of current model from groud truth
        %
        if optimValues.truth.dfval < options.TolFun
            if fitM && max(optimValues.truth.dM) < options.TolX
                if fitInit && optimValues.truth.dInitial < options.TolX
                    exitflag=1;
                    msg=['Success. trueloglike - loglike < ' num2str(options.TolFun)...
                        ' and distance from true model to fit model < ' num2str(options.TolX)];
                    return;
                else
                    exitflag=-4;
                    msg=['Not enough data? trueloglike - loglike < ' num2str(options.TolFun)...
                        ' and distance from true M to fit M > ' num2str(options.TolX)...
                        ' but not Initial'];
                    return;
                end
            elseif fitInit && optimValues.truth.dInitial < options.TolX
                exitflag=1;
                msg=['Success. trueloglike - loglike < ' num2str(options.TolFun)...
                    ' and distance from true model to fit model < ' num2str(options.TolX)];
                return;
            else
                exitflag=-3;
                msg=['Not enough data? trueloglike - loglike < ' num2str(options.TolFun)...
                    ' despite distance from true M to fit M > ' num2str(options.TolX)];
                return;
            end
        end
    end
    %
%decide whether or not to continue
%using difference of current model from previous one
    %
    if optimValues.fval > options.TolFun
        exitflag=1;
        msg=['Success. loglike > ' num2str(options.TolFun)];
        return;
    elseif max([fitInit*optimValues.prev.dInitial fitM*optimValues.prev.dM]) < options.TolX ...
            && abs(optimValues.prev.dfval) < options.TolFunChange
        if isnan(options.TolFun) && isempty(optimValues.truth)
            exitflag=1;
        else
            exitflag=-2;
        end
        msg=['Reached local maximum. Change in loglike < ' num2str(options.TolFunChange)...
            '. Size of change in model < ' num2str(options.TolX)];
        return;
    end
    %
%decide whether or not to continue
%using change in log likelihood on heldvack data
%i.e. test generalisation to prevent overfitting
    %
    if ~isempty(optimValues.holdback)
        if optimValues.holdback.dfval < options.MinLogLikeInc
            exitflag=-2;
            msg=['Overfitting. Change in loglike(held back data) < ' num2str(options.MinLogLikeInc)];
        end
    end

end

