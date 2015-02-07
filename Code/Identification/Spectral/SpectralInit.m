function [ modelobj ] = SpectralInit( simobj,w,varargin )
%SPECTRALINIT Summary of this function goes here
%   Detailed explanation goes here

[P1,P2,P3]=CalcObsProbs(simobj);

modelobj=SynapseIdModel.Rand(w,varargin{:});

Obs=MakeObs(modelobj);
Obsinv=(Obs'*Obs)\Obs';

Initial=P1*Obsinv;
% modelobj=modelobj.setInitial(Initial);


M=cell(1,simobj.NumPlast);

for a=1:length(P2)
    M{a}=zeros(length(P2{a}));
    for b=1:length(P3)
        M{a}=M{a}+P2{b}\squeeze(sum(P3{b,a},2));
    end
    M{a}=Obs*M{a}*Obsinv;
end

% modelobj=modelobj.setM(M);

% modelobj=modelobj + 0.05*SynapseIdModel.Rand(w);
modelobj=0.05*modelobj + SynapseIdModel('M',M,'Initial',Initial,'w',w);

modelobj=modelobj.Normalise;


end

