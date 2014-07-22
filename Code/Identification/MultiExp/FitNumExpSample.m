function [ num_exp,prob,ParticipationRatio,samples ] = FitNumExpSample( data,synapseOptions )
%[nnum_exp,P]=FITMULTIEXPSAMPLE(data,synapseOptions) number of exponentials nedded
%to fit distribution of data 
%   data           = cell of row vectors of dwell times. One row held back for
%       testing fits. Which row held back is rotated, then averaged over.
%   synapseOptions = SynapseOptimset
%   num_exp = number of exponentials needed
%   P       = matrix of participation ratios for coefficients of exponentials

maxnum=synapseOptions.MaxStates;

ParticipationRatio=zeros(synapseOptions.NumReps,synapseOptions.NumSample);

coefficient=zeros(maxnum,synapseOptions.NumSample);

for i=1:synapseOptions.NumReps
    
    tau=rand(1,maxnum);
    c=rand(1,maxnum);
    c=c/sum(c);
    init=[tau c(1:end-1)];

    %init = [ zeros(1,3)  1.0000    2.0000    3.0000   zeros(1,3)  0.2000    0.4000];
    
    samples = slicesample(init,synapseOptions.NumSample,...
        'logpdf',@JointDist);

    coefficient(1:maxnum-1,:)=samples(:,maxnum+1:end)';
    coefficient(maxnum,:)=1-sum(coefficient(1:maxnum-1,:),1);
    coefficient=abs(coefficient);
%     coefficient=sort(abs(coefficient),'descend');

    ParticipationRatio(i,:)=((sum(coefficient,1).^2)./sum(coefficient.^2,1));

end
% hist(P);
num_exp=round(mean(ParticipationRatio(:)));
count=hist(ParticipationRatio(:),1:synapseOptions.MaxStates);
prob=count(num_exp)/(numel(ParticipationRatio));


    function logjoint=JointDist(params)
        logjoint=-MultiExpLike(params',data)-MultiExpPrior(params',synapseOptions.PriorCcoeff,synapseOptions.PriorQcoeff);
    end

end

