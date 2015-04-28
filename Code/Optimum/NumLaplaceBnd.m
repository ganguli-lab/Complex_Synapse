function [ chains ] = NumLaplaceBnd( srange,nstates,trange,sym,varargin )
%chains=NUMLAPLACEBND(srange,nstates,trange,sym) numeric laplace bound
%   chains  = struct array (size=[1 length(srange)])
%   srange  = values of Laplace parameter at which we maximise
%   nstates = number of states in chain
%   trange  = values of time for snr curve
%   sym     = search for symmetric chains?
%   chains.s   = value of Laplace parameter at which we optimised
%   chains.qv  = nearest neighbour transitions of optimal model
%   chains.A   = value of Laplace transform at s for optimal model
%   chains.snr = snr curve of optimal model

chains(1,length(srange))=struct('s',[],'qv',[],'A',[],'snr',[],'KTp',[],'KTm',[]);

for i=1:length(srange)
    chains(i).s=srange(i);
    
    if sym
        [chains(i).qv,chains(i).A]=FindOptChainSL(srange(i),nstates,50,varargin{:});
        qp=chains(i).qv;
        qm=wrev(qp);
    else
        [chains(i).qv,chains(i).A]=FindOptChainL(srange(i),nstates,50,varargin{:});
        qp=chains(i).qv(1:nstates-1);
        qm=chains(i).qv(nstates:end);
    end
    
    [Wp,Wm,w]=MakeMultistate(qp,qm);
    modelobj=SynapseMemoryModel('Wp',Wp,'Wm',Wm,'w',w,'fp',0.5);
    
    chains(i).snr=modelobj.SNRcurve(trange);
    
    [~,dWp,dWm]=modelobj.SNRlaplaceGrad(srange(i));
    [chains(i).KTp,chains(i).KTm]=KTmults(Wp,Wm,dWp,dWm);
    
    DispCounter(i,length(srange),'');
end
DispCounter(length(srange)+1,length(srange),'');

end

