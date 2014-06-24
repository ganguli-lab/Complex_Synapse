function FitSynapseUtility_display( optimValues,options,state,exitflag )
%FITSYNAPSEUTILITY_UPDATE(optimValues,fitmodel,simobj,updaterfn,extraArgs)
%perform update of currrent model fit
%   optimValues = struct with information about the current state of the optimiser
%   options     = struct of options (see SynapseOptimset)
%   fitmodel    = SynapseIdModel where we'll store the model fit
%   simobj      = vector of SynapsePlastSeq with data to be fit
%   state       = string withh current state of updater
%   exitflag    = integer describing reason for exiting

switch options.Display
    case 'off'
        return;
    case 'iter'
        if strcmpi(state,'iter')
            DisplayMessage;
        end
    case 'final'
        if strcmpi(state,'done')
            DisplayMessage;
        end
    case 'notify'
        if strcmpi(state,'done') && exitflag < 1
            DisplayMessage;
        end
end
        
    function DisplayMessage
        disp(['Iteration: ' int2str(optimValues.iteration)]);
    end

end

