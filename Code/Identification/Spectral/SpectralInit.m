function [ modelobj ] = SpectralInit( seqobj,w,varargin )
%modelobj=SPECTRALINIT(seqobj,w,...) Crappy version of spectral algorithm
%for synapse identification. CAn be used as initialisation for EM.
%   modelobj = SynapseIdModel
%   seqobj = array of SynapsePlastSeq
%   w      = vector of synaptic weight labels for modelobj
%param/val:
%   RandFrac = coefficient of random model to add to modelobj (default=0.05)
%   others passed to SynapseIdModel.Rand

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='SpectralInit';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addParameter('RandFrac',0.05,@(x)validateattributes(x,{'numeric'},{'nonnegative'},'SpectralInit','options',3));
end
p.parse(varargin{:});


[P1,P2,P3]=CalcObsProbs(seqobj);

modelobj=SynapseIdModel.Rand(w,p.Unmatched);

Obs=MakeObs(modelobj);
Obsinv=(Obs'*Obs)\Obs';

Initial=P1*Obsinv;
% modelobj=modelobj.setInitial(Initial);


M=cell(1,seqobj.NumPlast);

for a=1:length(P2)
    M{a}=zeros(length(P2{a}));
    for b=1:length(P3)
        M{a}=M{a}+P2{b}\squeeze(sum(P3{b,a},2));
    end
    M{a}=Obs*M{a}*Obsinv;
end

% modelobj=modelobj.setM(M);

% modelobj=modelobj + 0.05*SynapseIdModel.Rand(w);
modelobj=p.Results.RandFrac*modelobj + SynapseIdModel('M',M,'Initial',Initial,'w',w);

modelobj=modelobj.Normalise;


end

