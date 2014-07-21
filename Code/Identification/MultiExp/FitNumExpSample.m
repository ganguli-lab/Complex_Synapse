function [ num_exp,prob,P ] = FitNumExpSample( data,synapseOptions )
%[nnum_exp,P]=FITMULTIEXPSAMPLE(data,synapseOptions) number of exponentials nedded
%to fit distribution of data 
%   data           = cell of row vectors of dwell times. One row held back for
%       testing fits. Which row held back is rotated, then averaged over.
%   synapseOptions = SynapseOptimset
%   num_exp = number of exponentials needed
%   P       = matrix of participation ratios for coefficients of exponentials

maxnum=synapseOptions.MaxStates;

P=zeros(synapseOptions.NumReps,synapseOptions.NumSample);

c=zeros(maxnum,synapseOptions.NumSample);

for i=1:synapseOptions.NumReps

    rnd = slicesample(rand(1,2*maxnum-1),synapseOptions.NumSample,...
        'logpdf',@JointDist);

    c(1:maxnum-1,:)=rnd(:,maxnum+1:end)';
    c(maxnum,:)=1-sum(c(1:maxnum-1,:),1);
    c=sort(abs(c),'descend');

    P(i,:)=((sum(c,1).^2)./sum(c.^2,1));

end
% hist(P);
num_exp=round(mean(P(:)));
count=hist(P(:),1:synapseOptions.MaxStates);
prob=count(num_exp)/(numel(P));


    function logjoint=JointDist(params)
        logjoint=-MultiExpLike(params',data)-MultiExpPrior(params',synapseOptions.PriorCcoeff,synapseOptions.PriorQcoeff);
    end

end

