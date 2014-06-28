function [ num_events,prob_st,KL,Ln ] = TryParams( truemodel,varargin )
%[num_events,prob_st,KL,Ln]=TRYPARAMS(truemodel,options,...) testing synapse
%fitting parameters
%   num_events = total number of plasticity events in simulations
%   prob_st    = probability of getting # states correct for this # events
%   KL         = mean and std of KL divergence from true model for this # events
%   Ln         = mean and std of Ln distance from true model for this # events
%   truemodel = SynapseIdModel used for simulations
%   options   = output of SynapseOptimset (optional)
%parameter/value pairs (with default)
%   num_t      = number of events per sequence (400)
%   min_ch     = minimum # sequences (10)
%   max_ch     = maximum # sequences (100)
%   crossval   = number of splits for ccross-validation in size fitting (2)
%   num_trials = number of repetitions for each # sequences (20)
%   plot       = plot results? (false)

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='TryParams';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addOptional('options',SynapseOptimset,@(x)validateattributes(x,{'struct'},{}));
    p.addParameter('num_t',400);
    p.addParameter('min_ch',10);
    p.addParameter('max_ch',100);
    p.addParameter('crossval',2);
    p.addParameter('num_trials',20);
    p.addParameter('plot',false);
%     p.addParameter('Normalise',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
end
p.parse(varargin{:});

fitsizeparams={'fp',truemodel.fp,'NormPower',truemodel.NormPower,'NormRows',truemodel.NormRows};

n_loop=floor(log2(p.Results.max_ch/p.Results.min_ch))+1;

num_events=zeros(1,n_loop);
prob_st=zeros(1,n_loop);
KL=zeros(2,n_loop);
Ln=zeros(2,n_loop);

num_ch=p.Results.min_ch;
loop=0;

while num_ch <= p.Results.max_ch
    fprintf('#seq: %d          ',num_ch);
    loop=loop+1;
    
    simobj=truemodel.Simulate(rand(2,p.Results.num_t,p.Results.crossval,num_ch));
    
    KLtrials=zeros(1,p.Results.num_trials);
    Lntrials=KLtrials;
    
    for trial=1:p.Results.num_trials
        fprintf([repmat('\b',[1 8+numel(int2str(trial-1))]) 'trial#: %d'] ,trial);
        
        fitmodel=FitSynapseSize(simobj,p.Results.options,fitsizeparams{:});
        prob_st(loop) = prob_st(loop) + truemodel.MatchW(fitmodel);
        
        guessmodel=truemodel.Randomise;
        fitmodel=FitSynapse(simobj(:)',guessmodel,p.Results.options);
        
        [~,KLfit]=truemodel.KLdivs(fitmodel);
        KLtrials(trial)=sum(KLfit);
        
        [~,Lnfit]=truemodel.LnNorm(fitmodel);
        Lntrials(trial)=sum(Lnfit);
        
    end%for trial
    
    num_events(loop)=p.Results.crossval*num_ch*p.Results.num_t;
   
    KL(1,loop)=mean(KLtrials);
    KL(2,loop)=std(KLtrials);
    
    Ln(1,loop)=mean(Lntrials);
    Ln(2,loop)=std(Lntrials);
        
    fprintf(repmat('\b',[1 15+numel(int2str(num_ch))+numel(int2str(trial))]))
    num_ch=2*num_ch;
end%while num_ch

prob_st=prob_st/p.Results.num_trials;

if p.Results.plot
    figure;
    subplot(1,3,1);
    plot(num_events,prob_st);
    ylim([0 1]);
    xlabel('# events');
    ylabel('prob # states correct');
    subplot(1,3,2);
    errorbar(num_events,KL(1,:),KL(2,:));
    xlabel('# events');
    ylabel('KL distance');
    subplot(1,3,3);
    errorbar(num_events,Ln(1,:),Ln(2,:));
    xlabel('# events');
    ylabel('L^n distance');
end

end

