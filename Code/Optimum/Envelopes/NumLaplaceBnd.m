function [ mats ] = NumLaplaceBnd( srange,nstates,mode,varargin )
%chains=NUMLAPLACEBND(srange,nstates,mode) numeric laplace bound
%   mats  = struct array (size=[1 length(srange)])
%   srange  = values of Laplace parameter at which we maximise
%   nstates = number of states in chain
%   trange  = values of time for snr curve
%   mode    = search for symmetric chains? include homeostatic plasticity?
%   mats.s   = value of Laplace parameter at which we optimised
%   mats.modelobj  = SynapseMemoryModel of optimal model
%   mats.A   = value of Laplace transform at s for optimal model
%   mats.snr = snr curve of optimal model

reps=50;

w=BinaryWeights(nstates);

mats(1,length(srange))=struct('s',[],'modelobj',[],'A',[],'snrb',[]);

for i=1:length(srange)
    
    DispCounter(i,length(srange),'s val: ');

    mats(i).s=srange(i);
    
    switch mode
        case 'hom'
            [Wp,Wm,Q,mats(i).A]=FindOptHomL(srange(i),nstates,reps,varargin{:});
            mats(i).modelobj=SynapseMemoryModel('Wp',Wp+Q,'Wm',Wm+Q,'w',w,'fp',0.5);
        otherwise
            [Wp,Wm,mats(i).A]=FindOptL(srange(i),nstates,reps,varargin{:});
            mats(i).modelobj=SynapseMemoryModel('Wp',Wp,'Wm',Wm,'w',w,'fp',0.5);
    end
    
    mats(i).snrb=mats(i).modelobj.SNRrunAve(1./srange);
end
DispCounter(length(srange)+1,length(srange),'s val: ');

end

