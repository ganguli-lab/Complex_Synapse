function [ env ] = NumLaplaceBndDouble( srange,nstates,sc,Ac,mode,varargin )
%chains=NUMLAPLACEBNDDOUBLE(srange,nstates,sc,Ac,mode) numeric laplace
%bound with constraint
%   env = struct
%   env.mats  = struct array (size=[1 length(srange)])
%   mats.s   = value of Laplace parameter at which we optimised
%   mats.modelobj  = SynapseMemoryModel of optimal model
%   mats.A   = value of Laplace transform at s for optimal model
%   mats.snr = snr curve of optimal model
%   srange  = values of Laplace parameter at which we maximise
%   nstates = number of states in chain
%   trange  = values of time for snr curve
%   sc,Ac   = constaint, A(sc)=Ac
%   mode    = search for symmetric chains? include homeostatic plasticity?


% function [ env ] = NumLaplaceBndDouble( srange,nstates,trange,sc,Ac,varargin )



if isnumeric(varargin{1})
    reps=varargin{1};
    varargin(1)=[];
else
    reps=50;
end

w=BinaryWeights(nstates);

% mats(1,length(srange))=struct('s',[],'modelobj',[],'A',[],'snr',[],'KTp',[],'KTm',[]);
mats(1,length(srange))=struct('s',[],'modelobj',[],'A',[],'snrb',[]);

for i=1:length(srange)
    
    DispCounter(i,length(srange),'s val: ');

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
%     
%     [~,dWp,dWm]=mats(i).modelobj.SNRlaplaceGrad(srange(i));
%     [mats(i).KTp,mats(i).KTm]=KTmults(mats(i).modelobj.Wp,mats(i).modelobj.Wm,dWp,dWm);
end
DispCounter(length(srange)+1,length(srange),'s val: ');


Aenv=[mats.A];
Aenv=Aenv.*srange;

env=struct('sc',sc,'Ac',Ac,'mats',mats,'tau',1./srange,'SNRbenv',Aenv);

end

