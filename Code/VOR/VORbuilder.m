function [ exptobj ] = VORbuilder( builder_h,numStates,paramPot,paramWT,paramKO,fpNorm,fpInc,fpDec,t_inc,t_dec )
%VORBUILDER Summary of this function goes here
%   Detailed explanation goes here



[Wp,~,w]=builder_h(numStates,paramPot);
[~,WmWT]=builder_h(numStates,paramWT);
[~,WmKO]=builder_h(numStates,paramKO);

WTobj=SynapseMemoryModel('Wp',Wp,'Wm',WmWT,'w',w,'fp',fpNorm);
KOobj=SynapseMemoryModel('Wp',Wp,'Wm',WmKO,'w',w,'fp',fpNorm);

nopreobj=VORtrainSeq('tTrain',t_dec+t_inc,'fps',[fpNorm fpInc]);
withpreobj=VORtrainSeq('tTrain',[t_dec t_dec+t_inc],'fps',[fpNorm fpDec fpInc]);

exptobj=VORexperiment('WT',WTobj,'KO',KOobj,'nopre',nopreobj,'withpre',withpreobj);

end

