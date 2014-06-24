function [ optimValues ] = FitSynapseUtility_calcchanges( optimValues,options )
%[optimValues]=FITSYNAPSEUTILITY_CALCCHANGES(optimValues,options)
%calculate changes in model fit due to last update 
%   optimValues = struct with information about the current state of the optimiser
%   options     = struct of options (see SynapseOptimset)


%calculate size of changes in model
    %
    switch options.ModelDiff
        case 'KL'
            [optimValues.prev.dInitial,optimValues.prev.dM]=optimValues.prev.model.KLdivs(optimValues.model);
        case 'Ln'
            [optimValues.prev.dInitial,optimValues.prev.dM]=optimValues.prev.model.LnNorm(options.NormPower,optimValues.model);
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
                [optimValues.truth.dInitial,optimValues.truth.dM]=optimValues.truth.model.KLdivs(optimValues.model);
            case 'Ln'
                [optimValues.truth.dInitial,optimValues.truth.dM]=optimValues.truth.model.LnNorm(options.NormPower,optimValues.model);
        end
        %
    end

end

