function [ optimValues ] = FitSynapseUtility_calcchanges( optimValues,fitmodel,options )
%[optimValues]=FITSYNAPSEUTILITY_CALCCHANGES(optimValues,fitmodel,options)
%calculate changes in model fit due to last update 
%   optimValues = struct with information about the current state of the optimiser
%   fitmodel    = SynapseIdModel where we'll store the model fit
%   options     = struct of options (see SynapseOptimset)


%calculate size of changes in model
    %
    switch options.ModelDiff
        case 'KL'
            [optimValues.prev.dInitial,optimValues.prev.dM]=optimValues.prev.model.KLdivs(fitmodel);
        case 'Ln'
            [optimValues.prev.dInitial,optimValues.prev.dM]=optimValues.prev.model.LnNorm(options.NormPower,fitmodel);
    end
    optimValues.prev.dfval=optimValues.fval-optimValues.prev.fval;
    optimValues.funcCount=optimValues.funcCount+1;
    %
%calculate difference of current model from groud truth
    %
    if ~isempty(optimValues.truth)
        %
        optimValues.truth.dfval=optimValues.truth.fval-optimValues.fval;
        switch options.ModelDiff
            case 'KL'
                [optimValues.truth.dInitial,optimValues.truth.dM]=optimValues.truth.model.KLdivs(fitmodel);
            case 'Ln'
                [optimValues.truth.dInitial,optimValues.truth.dM]=optimValues.truth.model.LnNorm(options.NormPower,fitmodel);
        end
        %
    end

end

