function [ mats ] = NumLaplaceBnd( srange,nstates,trange,mode,varargin )
%chains=NUMLAPLACEBND(srange,nstates,trange,mode) numeric laplace bound
%   mats  = struct array (size=[1 length(srange)])
%   srange  = values of Laplace parameter at which we maximise
%   nstates = number of states in chain
%   trange  = values of time for snr curve
%   mode    = search for symmetric chains? include homeostatic plasticity?
%   chains.s   = value of Laplace parameter at which we optimised
%   chains.qv  = nearest neighbour transitions of optimal model
%   chains.A   = value of Laplace transform at s for optimal model
%   chains.snr = snr curve of optimal model

reps=50;

mats(1,length(srange))=struct('s',[],'Wp',[],'Wm',[],'Q',[],'A',[],'snr',[],'KTp',[],'KTm',[]);

for i=1:length(srange)
    
    DispCounter(i,length(srange),'s val: ');

    mats(i).s=srange(i);
    
    switch mode
        case 'hom'
            [mats(i).Wp,mats(i).Wm,mats(i).Q,mats(i).A]=FindOptHomL(srange(i),nstates,reps,varargin{:});
            Wp=mats(i).Wp+mats(i).Q;
            Wm=mats(i).Wm+mats(i).Q;
        otherwise
            [mats(i).Wp,mats(i).Wm,mats(i).A]=FindOptL(srange(i),nstates,reps,varargin{:});
            Wp=mats(i).Wp;
            Wm=mats(i).Wm;
    end
    
    w=BinaryWeights(nstates);
    modelobj=SynapseMemoryModel('Wp',Wp,'Wm',Wm,'w',w,'fp',0.5);
    
    mats(i).snr=modelobj.SNRcurve(trange);
    
    [~,dWp,dWm]=modelobj.SNRlaplaceGrad(srange(i));
    [mats(i).KTp,mats(i).KTm]=KTmults(modelobj.Wp,modelobj.Wm,dWp,dWm);
end
DispCounter(length(srange)+1,length(srange),'s val: ');

end

