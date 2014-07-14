function [ num_exp,P ] = FitNumExpSample( data,synapseOptions )
%[nnum_exp,P]=FITMULTIEXPSAMPLE(data,synapseOptions) number of exponentials nedded
%to fit distribution of data 
%   data           = cell of row vectors of dwell times. One row held back for
%       testing fits. Which row held back is rotated, then averaged over.
%   synapseOptions = SynapseOptimset
%   num_exp = number of exponentials needed
%   P       = matrix of participation ratios for coefficients of exponentials

P=zeros(synapseOptions.NumReps,synapseOptions.NumSample);

for i=1:synapseOptions.NumReps

    maxnum=synapseOptions.MaxStates;
    rnd = slicesample(rand(1,2*maxnum-1),synapseOptions.NumSample,'logpdf',...
        @JointDist);

    c=rnd(:,maxnum+1:end)';
    c=[c;1-sum(c,1)];
    c=sort(abs(c),'descend');

    P(i,:)=((sum(c,1).^2)./sum(c.^2,1));

end
% hist(P);
num_exp=round(mean(P(:)));


    function logjoint=JointDist(params)
        logjoint=-MultiExpFun(params',data)-MultiExpPrior(params',synapseOptions.PriorCcoeff,synapseOptions.PriorQcoeff);
    end

end

