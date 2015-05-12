function [ mats ] = NumLaplaceBnd( srange,nstates,trange,varargin )
%chains=NUMLAPLACEBND(srange,nstates,trange,sym) numeric laplace bound
%   mats  = struct array (size=[1 length(srange)])
%   srange  = values of Laplace parameter at which we maximise
%   nstates = number of states in chain
%   trange  = values of time for snr curve
%   sym     = search for symmetric chains?
%   hom     = include homeostatic plasticity?
%   chains.s   = value of Laplace parameter at which we optimised
%   chains.qv  = nearest neighbour transitions of optimal model
%   chains.A   = value of Laplace transform at s for optimal model
%   chains.snr = snr curve of optimal model

mats(1,length(srange))=struct('s',[],'Wp',[],'Wm',[],'Q',[],'A',[],'snr',[],'KTp',[],'KTm',[]);

for i=1:length(srange)
    
    DispCounter(i,length(srange),'s val: ');

    mats(i).s=srange(i);
    
    [mats(i).Wp,mats(i).Wm,mats(i).Q,mats(i).A]=FindOptHomL(srange(i),nstates,50,varargin{:});
    
    w=BinaryWeights(nstates);
    modelobj=SynapseMemoryModel('Wp',mats(i).Wp+mats(i).Q,'Wm',mats(i).Wm+mats(i).Q,'w',w,'fp',0.5);
    
    mats(i).snr=modelobj.SNRcurve(trange);
    
    [~,dWp,dWm]=modelobj.SNRlaplaceGrad(srange(i));
    [mats(i).KTp,mats(i).KTm]=KTmults(modelobj.Wp,modelobj.Wm,dWp,dWm);
end
DispCounter(length(srange)+1,length(srange),'s val: ');

end

