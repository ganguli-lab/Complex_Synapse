function [ chains ] = FixNumLaplaceBndChain( oldchains,srange,nstates,trange,inds,mode,varargin )
%chains=FIXNUMLAPLACEBNDCHAIN(srange,nstates,trange,sym) numeric laplace bound
%   chains  = struct array (size=[1 length(srange)])
%   srange  = values of Laplace parameter at which we maximise
%   nstates = number of states in chain
%   trange  = values of time for snr curve
%   sym     = search for symmetric chains?
%   hom     = include homeostatic plasticity?
%   chains.s   = value of Laplace parameter at which we optimised
%   chains.qv  = nearest neighbour transitions of optimal model
%   chains.A   = value of Laplace transform at s for optimal model
%   chains.snr = snr curve of optimal model

chains=oldchains;
reps=200;

for j=1:length(inds)
    i=inds(j);
    
    DispCounter(j,length(inds),'s val: ');
    chains(i).s=srange(i);
    
    switch mode
        case 'n'
            [chains(i).qv,chains(i).A]=FindOptChainL(srange(i),nstates,reps,varargin{:});
            qp=chains(i).qv(1:nstates-1);
            qm=chains(i).qv(nstates:end);
        case 's'
            [chains(i).qv,chains(i).A]=FindOptChainHomL(srange(i),nstates,reps,varargin{:});
            [qp,qm]=MakeHomq(chains(i).qv,0.5);
        case 'c'
            [chains(i).qv,chains(i).A]=FindOptChainHomLC(srange(i),nstates,reps,varargin{:});
            [qp,qm]=MakeHomqC(chains(i).qv,0.5);
        case 'a'
            [chains(i).qv,chains(i).A]=FindOptChainHomLA(srange(i),nstates,reps,varargin{:});
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

