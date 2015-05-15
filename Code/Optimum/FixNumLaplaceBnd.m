function [ mats ] = FixNumLaplaceBnd( oldmats,srange,nstates,trange,inds,varargin )
%chains=FIXNUMLAPLACEBND(srange,nstates,trange,sym) numeric laplace bound
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

mats=oldmats;
reps=200;

for j=1:length(inds)
    i=inds(j);
    
    DispCounter(j,length(inds),'s val: ');
    mats(i).s=srange(i);
    
    [mats(i).Wp,mats(i).Wm,mats(i).Q,mats(i).A]=FindOptHomL(srange(i),nstates,reps,varargin{:});
    
    w=BinaryWeights(nstates);
    modelobj=SynapseMemoryModel('Wp',mats(i).Wp+mats(i).Q,'Wm',mats(i).Wm+mats(i).Q,'w',w,'fp',0.5);
    
    mats(i).snr=modelobj.SNRcurve(trange);
    
    [~,dWp,dWm]=modelobj.SNRlaplaceGrad(srange(i));
    [mats(i).KTp,mats(i).KTm]=KTmults(modelobj.Wp,modelobj.Wm,dWp,dWm);
    
    if mats(i).A < oldmats(i).A
        mats(i)=oldmats(i);
    end
end
DispCounter(length(inds)+1,length(inds),'s val: ');

end

