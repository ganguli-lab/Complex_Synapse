function [ num_exp,fvals ] = FitNumExp( data,synapseOptions,optimOptions )
%FITNUMEXP Summary of this function goes here
%   Detailed explanation goes here


num_exp=0;
negloglike=Inf;
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');

while num_exp < synapseOptions.MaxStates
    
    fitnegloglike=zeros(1,length(data));
    holdinds=1:length(data);
    for holdback=holdinds
        fitinds=holdinds;
        fitinds(holdback)=[];
        fitnegloglike(holdback)=OneCrossVal(num_exp+1,[data{fitinds}],data{holdback});
    end%for holdback
    fitnegloglike=mean(fitnegloglike);
    
    if negloglike-fitnegloglike > synapseOptions.MinLogLikeInc
        num_exp=num_exp+1;
        negloglike=fitnegloglike;
        fvals(num_exp)=negloglike;
    else
        break;
    end%if fit improves enough
    
end%while num_exp

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

end

