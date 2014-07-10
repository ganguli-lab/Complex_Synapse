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

tic;

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='TryParams';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addOptional('synapseOptions',SynapseOptimset,@(x)validateattributes(x,{'SynapseOptimset'},{},'TryParams','synapseOptions',2));
    p.addOptional('optimOptions',optimoptions('fmincon','Algorithm','interior-point'),@(x)validateattributes(x,{'optim.options.Fmincon'},{},'TryParams','optimOptions',3));
    p.addParameter('num_t',400);
    p.addParameter('min_ch',10);
    p.addParameter('max_ch',100);
    p.addParameter('crossval',2);
    p.addParameter('num_trials',20);
    p.addParameter('plot',false);
    p.addParameter('doSize',true);
    p.addParameter('doDist',true);
%     p.addParameter('Normalise',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
end
p.parse(varargin{:});
r=p.Results;

fitsizeparams={'fp',truemodel.fp,'NormPower',truemodel.NormPower,'NormRows',truemodel.NormRows};

n_loop=floor(log2(p.Results.max_ch/p.Results.min_ch))+1;

num_events=zeros(1,n_loop);
prob_st=zeros(2,n_loop);
KL=zeros(2,n_loop);
Ln=zeros(2,n_loop);

num_ch=round(r.min_ch/r.crossval);
loop=0;

while num_ch*r.crossval <= r.max_ch
    fprintf('#seq: %d          ',num_ch);
    loop=loop+1;
    
    simobj=truemodel.Simulate(rand(2,r.num_t,r.crossval,num_ch));
    
    KLtrials=zeros(1,r.num_trials);
    Lntrials=KLtrials;
    
    for trial=1:r.num_trials
        fprintf([repmat('\b',[1 8+numel(int2str(trial-1))]) 'trial#: %d'] ,trial);
        
        if r.doSize
        
%             fitmodel=FitSynapseSize(simobj,r.synapseOptions,fitsizeparams);
            fitmodel=FitSynapseSizeDwell(simobj,r.synapseOptions,r.optimOptions,fitsizeparams);
            prob_st(1,loop) = prob_st(loop) + truemodel.MatchW(fitmodel);
        
        end
        
        if r.doDist
            
            guessmodel=truemodel.Randomise;
            fitmodel=FitSynapse(simobj(:)',guessmodel,r.synapseOptions);

            [~,KLfit]=truemodel.KLdivs(fitmodel);
            KLtrials(trial)=sum(KLfit);

            [~,Lnfit]=truemodel.LnNorm(fitmodel);
            Lntrials(trial)=sum(Lnfit);
        
        end
        
    end%for trial
    
    num_events(loop)=r.crossval*num_ch*r.num_t;
   
    KL(1,loop)=mean(KLtrials);
    KL(2,loop)=std(KLtrials);
    
    Ln(1,loop)=mean(Lntrials);
    Ln(2,loop)=std(Lntrials);
        
    fprintf(repmat('\b',[1 15+numel(int2str(num_ch))+numel(int2str(trial))]))
    num_ch=2*num_ch;
end%while num_ch

prob_st=prob_st/r.num_trials;
prob_st(2,:) = sqrt( prob_st(1,:) .* (1-prob_st(1,:)) / r.num_trials );

toc;

if r.plot
    figure;
    if r.doSize
        subplot(1,1+2*r.doDist,1);
        errorbar(num_events,prob_st(1,:),prob_st(2,:));
        ylim([0 1]);
        xlabel('# events');
        ylabel('prob # states correct');
    end
    if r.doDist
        subplot(1,r.doSize+2,r.doSize+1);
        errorbar(num_events,KL(1,:),KL(2,:));
        xlabel('# events');
        ylabel('KL distance');
        subplot(1,r.doSize+2,r.doSize+2);
        errorbar(num_events,Ln(1,:),Ln(2,:));
        xlabel('# events');
        ylabel('L^n distance');
    end
end

end

