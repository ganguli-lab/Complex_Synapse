function [ chains ] = NumLaplaceBndChain( srange,nstates,mode,varargin )
%chains=NUMLAPLACEBNDCHAIN(srange,nstates,mode) numeric laplace bound
%   chains  = struct array (size=[1 length(srange)])
%   srange  = values of Laplace parameter at which we maximise
%   nstates = number of states in chain
%   mode    = search for symmetric chains? include homeostatic plasticity?
%   chains.s   = value of Laplace parameter at which we optimised
%   chains.qv  = nearest neighbour transitions of optimal model
%   chains.A   = value of Laplace transform at s for optimal model
%   chains.snr = snr curve of optimal model

chains(1,length(srange))=struct('s',[],'qv',[],'A',[],'modelobj',[],'snrb',[]);
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
    chains(i).modelobj=SynapseMemoryModel('Wp',Wp,'Wm',Wm,'w',w,'fp',0.5);
    
    chains(i).snrb=chains(i).modelobj.SNRrunAve(1./srange);
    
end%for i
DispCounter(length(srange)+1,length(srange),'s val: ');

end

