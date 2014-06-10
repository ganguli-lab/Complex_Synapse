function [ optimValues ] = FitSynapseUtility_calcchanges( optimValues,fitmodel )
%[optimValues]=FITSYNAPSEUTILITY_CALCCHANGES(optimValues,fitmodel) calculate changes
%in model fit due to last update
%   optimValues = struct with information about the current state of the optimiser
%   fitmodel    = SynapseIdModel where we'll store the model fit


%calculate size of changes in model
    %
    [optimValues.prev.KLInitial,optimValues.prev.KLM]=optimValues.prev.model.KLdivs(fitmodel);
    optimValues.prev.dfval=optimValues.fval-optimValues.prev.fval;
    optimValues.funcCount=optimValues.funcCount+1;
    %
%calculate difference of current model from groud truth
    %
    if ~isempty(optimValues.truth)
        %
        optimValues.truth.dfval=optimValues.truth.fval-optimValues.fval;
        [optimValues.truth.KLInitial,optimValues.truth.KLM]=optimValues.truth.model.KLdivs(fitmodel);
        %
    end

end

