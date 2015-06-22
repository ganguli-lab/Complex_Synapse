function [ chains ] = FixNumLaplaceBndChain( oldchains,srange,nstates,trange,inds,mode,varargin )
%chains=FIXNUMLAPLACEBNDCHAIN(srange,nstates,trange,inds,mode) numeric laplace bound
%   chains  = struct array (size=[1 length(srange)])
%   srange  = values of Laplace parameter at which we maximise
%   nstates = number of states in chain
%   trange  = values of time for snr curve
%   mode    = search for symmetric chains? include homeostatic plasticity?
%   chains.s   = value of Laplace parameter at which we optimised
%   chains.qv  = nearest neighbour transitions of optimal model
%   chains.A   = value of Laplace transform at s for optimal model
%   chains.snr = snr curve of optimal model

chains=oldchains;
reps=1000;

for j=1:length(inds)
    i=inds(j);
    
    DispCounter(j,length(inds),'s val: ');
    chains(i).s=srange(i);
    
    switch mode
        case 'sym'
            [chains(i).qv,chains(i).A]=FindOptChainSL(srange(i),nstates,50,varargin{:});
            qp=chains(i).qv;
            qm=wrev(qp);
        case 'homs'
            [chains(i).qv,chains(i).A]=FindOptChainHomL(srange(i),nstates,reps,varargin{:});
            [qp,qm]=MakeHomq(chains(i).qv,0.5);
        case 'homc'
            [chains(i).qv,chains(i).A]=FindOptChainHomLC(srange(i),nstates,reps,varargin{:});
            [qp,qm]=MakeHomqC(chains(i).qv,0.5);
        case 'homa'
            [chains(i).qv,chains(i).A]=FindOptChainHomLA(srange(i),nstates,reps,varargin{:});
            qp=chains(i).qv(1:nstates-1);
            qm=chains(i).qv(nstates:end);
        otherwise
            [chains(i).qv,chains(i).A]=FindOptChainL(srange(i),nstates,reps,varargin{:});
            qp=chains(i).qv(1:nstates-1);
            qm=chains(i).qv(nstates:end);
    end
    
    [Wp,Wm,w]=MakeMultistate(qp,qm);
    modelobj=SynapseMemoryModel('Wp',Wp,'Wm',Wm,'w',w,'fp',0.5);
    
    chains(i).snr=modelobj.SNRcurve(trange);
    
    [~,dWp,dWm]=modelobj.SNRlaplaceGrad(srange(i));
    [chains(i).KTp,chains(i).KTm]=KTmults(Wp,Wm,dWp,dWm);
    
    if chains(i).A < oldchains(i).A
        chains(i)=oldchains(i);
    end
end
DispCounter(length(inds)+1,length(inds),'s val: ');

end

