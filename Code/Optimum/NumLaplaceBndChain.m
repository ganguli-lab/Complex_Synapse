function [ chains ] = NumLaplaceBndChain( srange,nstates,trange,mode,varargin )
%chains=NUMLAPLACEBNDCHAIN(srange,nstates,trange,sym) numeric laplace bound
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

chains(1,length(srange))=struct('s',[],'qv',[],'A',[],'snr',[],'KTp',[],'KTm',[]);
reps=50;

for i=1:length(srange)
    
    DispCounter(i,length(srange),'s val: ');

    chains(i).s=srange(i);
    
    switch mode
        case 'sym'
            [chains(i).qv,chains(i).A]=FindOptChainSL(srange(i),nstates,reps,varargin{:});
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
    end%switch mode
    
    [Wp,Wm,w]=MakeMultistate(qp,qm);
    modelobj=SynapseMemoryModel('Wp',Wp,'Wm',Wm,'w',w,'fp',0.5);
    
    chains(i).snr=modelobj.SNRcurve(trange);
    
    [~,dWp,dWm]=modelobj.SNRlaplaceGrad(srange(i));
    [chains(i).KTp,chains(i).KTm]=KTmults(Wp,Wm,dWp,dWm);
end%for i
DispCounter(length(srange)+1,length(srange),'s val: ');

end

