function [ env ] = FixNumLaplaceBndDouble( oldenv,srange,nstates,sc,Ac,mode,inds,varargin )
%env=FIXNUMLAPLACEBNDDOUBLE(oldenv,srange,nstates,sc,Ac,,mode,inds) numeric laplace bound
%   env.mats    = struct array (size=[1 length(srange)])
%   srange  = values of Laplace parameter at which we maximise
%   nstates = number of states in chain
%   sc,Ac   = constaint, A(sc)=Ac
%   mode    = search for symmetric chains? include homeostatic plasticity?
%   mats.s  = value of Laplace parameter at which we optimised
%   mats.modelobj = SynapseMemoryModel of optimal model
%   mats.A    = value of Laplace transform at s for optimal model
%   mats.snrb = snr curve of optimal model

mats=oldenv.mats;
if isnumeric(varargin{1})
    reps=varargin{1};
    varargin(1)=[];
else
    reps=200;
end
% reps=50;
w=BinaryWeights(nstates);

for j=1:length(inds)
    i=inds(j);
    
    DispCounter(j,length(inds),'s val: ');
    mats(i).s=srange(i);
    
    switch mode
        case 'ineq'
            [Wp,Wm,mats(i).A]=FindOptDoubleIL(srange(i),sc,Ac,nstates,reps,varargin{:});
            mats(i).modelobj=SynapseMemoryModel('Wp',Wp,'Wm',Wm,'w',w,'fp',0.5);
        otherwise
            [Wp,Wm,mats(i).A]=FindOptDoubleL(srange(i),sc,Ac,nstates,reps,varargin{:});
            mats(i).modelobj=SynapseMemoryModel('Wp',Wp,'Wm',Wm,'w',w,'fp',0.5);
    end
    
    if mats(i).modelobj.isvalid
        mats(i).snrb=mats(i).modelobj.SNRrunAve(1./srange);
    end
    
    if mats(i).A < oldenv.mats(i).A
        mats(i)=oldenv.mats(i);
    end
end
DispCounter(length(inds)+1,length(inds),'s val: ');


Aenv=[mats.A];
Aenv=Aenv.*srange;

env=struct('sc',sc,'Ac',Ac,'mats',mats,'tau',1./srange,'SNRbenv',Aenv);


end

