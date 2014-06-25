function FitSynapseUtility_display( optimValues,options,state,exitflag,msg )
%FITSYNAPSEUTILITY_UPDATE(optimValues,fitmodel,simobj,updaterfn,extraArgs)
%perform update of currrent model fit
%   optimValues = struct with information about the current state of the optimiser
%   options     = struct of options (see SynapseOptimset)
%   fitmodel    = SynapseIdModel where we'll store the model fit
%   simobj      = vector of SynapsePlastSeq with data to be fit
%   state       = string withh current state of updater
%   exitflag    = integer describing reason for exiting
%   msg         = string describing reason for exiting

switch options.Display
    case 'off'
        return;
    case 'addstate'
        return;
    case 'iter'
        if strcmpi(state,'iter')
            DisplayMessage;
            DisplayMessageDiff;
            DisplayMessageDiffTruth;
        end
        if strcmpi(state,'init')
            DisplayMessageTruth;
        end
        if strcmpi(state,'done')
            DisplayMessage;
            disp(msg);
        end
    case 'final'
        if strcmpi(state,'done')
            disp(msg);
            DisplayMessage;
        end
    case 'notify'
        if strcmpi(state,'done') && exitflag < 1
            disp(msg);
            DisplayMessage;
            DisplayMessageTruth;
            DisplayMessageDiffTruth;
        end
end
        
    function DisplayMessage
        disp(['Iterations: ' num2str(optimValues.iteration,3) '  LogLike: ' num2str(optimValues.fval,'%6.2f')]);
    end

    function DisplayMessageDiff
      disp([' Modeldiff: LogLike: ' num2str(optimValues.prev.dfval,'%6.2f') '  Initial: ' num2str(optimValues.prev.dInitial,'%4.1e') '  M: ' num2str(optimValues.prev.dM,'% 4.1e')]);
    end

    function DisplayMessageDiffTruth
        if ~isempty(optimValues.truth)
            disp([' True Modeldiff: LogLike: '  num2str(optimValues.truth.dfval,'%6.2f')  ' Initial: ' num2str(optimValues.truth.dInitial,'%4.1e') ' M: ' num2str(optimValues.truth.dM,'% 4.1e')]);
        end
    end

    function DisplayMessageTruth
        if ~isempty(optimValues.truth)
            disp(['True LogLike: ' num2str(optimValues.truth.fval,'%6.2f')]);
        end
    end

end

