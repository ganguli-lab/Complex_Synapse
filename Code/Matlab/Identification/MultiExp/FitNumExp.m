function [ num_exp,fvals ] = FitNumExp( data,synapseOptions,optimOptions )
%[num_exp,fvals]=FITNUMEXP(data,synapseOptions,optimOptions) number of
%exponentials nedded to fit distribution of data
%   data           = cell of row vectors of dwell times. One row held back for
%       testing fits. Which row held back is rotated, then averaged over.
%   synapseOptions = SynapseOptimset
%   optimOptions   = optimoptions for fminunc
%   num_exp = number of exponentials needed
%   fvals   = row of neg-log-likelihoods for each value of num_exp tested


num_exp=0;
negloglike=Inf;
fvals=NaN(1,synapseOptions.MaxStates);
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');

while num_exp <= synapseOptions.MaxStates
    
    fitnegloglike=zeros(1,length(data));
    holdinds=1:length(data);
    for holdback=holdinds
        fitinds=holdinds;
        fitinds(holdback)=[];
        fitnegloglike(holdback)=OneCrossVal(num_exp+1,[data{fitinds}],data{holdback});
    end%for holdback
    fitnegloglike=mean(fitnegloglike);
    
    fvals(num_exp+1)=fitnegloglike;
    if negloglike-fitnegloglike > synapseOptions.MinLogLikeInc
        num_exp=num_exp+1;
        negloglike=fitnegloglike;
        DispStates;
    else
        break;
    end%if fit improves enough
    
end%while num_exp

num_exp=max(1,num_exp);

warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:singularMatrix');


    function crossnegloglike=OneCrossVal(trialnum_exp,fitdata,testdata)
        bestfit=rand(2*trialnum_exp-11,1);
        bestnegloglike=Inf;
        for rep=1:synapseOptions.NumReps
            try
                [x,xloglike]=FitMultiExp(trialnum_exp,fitdata,optimOptions);
                if xloglike<bestnegloglike
                    bestnegloglike=xloglike;
                    bestfit=x;
                end%if better fit
            catch ME
            end
        end%for i
        crossnegloglike=MultiExpFun(bestfit,testdata);
    end

    function DispStates
        if ~strcmpi(synapseOptions.Display,'off')
            disp(['# states: ' int2str(num_exp)]);
        end
    end



end

