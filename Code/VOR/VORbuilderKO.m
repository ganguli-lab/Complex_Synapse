function [ exptobj ] = VORbuilderKO( builder_h,numStates,paramPot,paramWT,paramKO,fpNorm,fpInc,fpDec,t_inc,t_dec,pooled )
%VORexptobj=VORBUILDERKO(builder_h,numStates,paramPot,paramWT,paramKO,fpNorm,fpInc,fpDec,t_inc,t_dec,pooled)
%building a VORexptKO object
%   builder_h  = function handle [Wp,Wm,w]=builder_h(numStates,param)
%   num_states = number of synaptic states
%   paramPot   = paramter for Wp for both WT/KO
%   paramWT    = paramter for Wm for WT
%   paramKO    = paramter for Wm for KO
%   fpNorm     = fraction of potentiating events at baseline
%   fpInc      = fraction of potentiating events during gain-increase training
%   fpDec      = fraction of potentiating events during gain-decrease pre-training
%   t_inc      = duration of gain-increase training
%   t_dec      = duration of gain-decrease pre-training
%   pooled     = is it a pooled resource model?



[Wp,~,w]=builder_h(numStates,paramPot);
[~,WmWT]=builder_h(numStates,paramWT);
[~,WmKO]=builder_h(numStates,paramKO);

WTobj=SynapseMemoryModel('Wp',Wp,'Wm',WmWT,'w',w,'fp',fpNorm);
KOobj=SynapseMemoryModel('Wp',Wp,'Wm',WmKO,'w',w,'fp',fpNorm);

nopreobj=VORtrainSeq('tTrain',t_inc,'fps',[fpNorm fpInc]);
withpreobj=VORtrainSeq('tTrain',[t_dec t_dec+t_inc],'fps',[fpNorm fpDec fpInc]);

exptobj=VORexptKO('WT',WTobj,'KO',KOobj,'nopre',nopreobj,'withpre',withpreobj,'pooled',pooled);

end

