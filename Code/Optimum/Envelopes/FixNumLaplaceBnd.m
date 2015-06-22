function [ mats ] = FixNumLaplaceBnd( oldmats,srange,nstates,trange,mode,inds,varargin )
%chains=FIXNUMLAPLACEBND(srange,nstates,trange,mode,inds) numeric laplace bound
%   chains  = struct array (size=[1 length(srange)])
%   srange  = values of Laplace parameter at which we maximise
%   nstates = number of states in chain
%   trange  = values of time for snr curve
%   mode    = search for symmetric chains? include homeostatic plasticity?
%   chains.s   = value of Laplace parameter at which we optimised
%   chains.qv  = nearest neighbour transitions of optimal model
%   chains.A   = value of Laplace transform at s for optimal model
%   chains.snr = snr curve of optimal model

mats=oldmats;
reps=200;
w=BinaryWeights(nstates);

for j=1:length(inds)
    i=inds(j);
    
    DispCounter(j,length(inds),'s val: ');
    mats(i).s=srange(i);
    
    switch mode
        case 'hom'
            [Wp,Wm,Q,mats(i).A]=FindOptHomL(srange(i),nstates,reps,varargin{:});
            mats(i).modelobj=SynapseMemoryModel('Wp',Wp+Q,'Wm',Wm+Q,'w',w,'fp',0.5);
        otherwise
            [Wp,Wm,mats(i).A]=FindOptL(srange(i),nstates,reps,varargin{:});
            mats(i).modelobj=SynapseMemoryModel('Wp',Wp,'Wm',Wm,'w',w,'fp',0.5);
    end
    
    mats(i).snr=mats(i).modelobj.SNRcurve(trange);
    
    [~,dWp,dWm]=mats(i).modelobj.SNRlaplaceGrad(srange(i));
    [mats(i).KTp,mats(i).KTm]=KTmults(mats(i).modelobj.Wp,mats(i).modelobj.Wm,dWp,dWm);
    
    if mats(i).A < oldmats(i).A
        mats(i)=oldmats(i);
    end
end
DispCounter(length(inds)+1,length(inds),'s val: ');

end

