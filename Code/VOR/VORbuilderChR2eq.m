function [ exptobj ] = VORbuilderChR2eq( builder_h,numStates,paramPot,paramDep,fpNorm,fpInc,fpNS,t_inc,pooled )
%VORexptobj=VORBUILDERCHR2EQ(builder_h,numStates,paramPot,paramDep,fpNorm,fpInc,fpNS,t_inc,pooled)
%building a VORexperiment object with non-specific ChR2 CF stim, assuming
%it is long enough to equilibriate
%   builder_h  = function handle [Wp,Wm,w]=builder_h(numStates,param)
%   num_states = number of synaptic states
%   paramPot   = paramter for Wp
%   paramDep   = paramter for Wm
%   fpNorm     = fraction of potentiating events at baseline
%   fpInc      = fraction of potentiating events during gain-increase training
%   fpNS       = fraction of potentiating events during non-specific CF stim
%   t_inc      = duration of gain-increase training
%   pooled     = is it a pooled resource model?



[Wp,~,w] = builder_h(numStates,paramPot);
[~,Wm] = builder_h(numStates,paramDep);

WTobj = SynapseMemoryModel('Wp',Wp,'Wm',Wm,'w',w,'fp',fpNorm);

VORrel = VORtrainSeq('tTrain',t_inc, 'fps',[fpNorm fpInc]);
VORcomp = VORtrainSeq('tTrain',t_inc, 'fps',[fpNorm fpNorm]);
nopreobj = VORtrainSeqDiff('VORrel',VORrel, 'VORcomp',VORcomp);

VORrel = VORtrainSeq('tTrain',t_inc, 'fps',[fpNS fpInc]);
VORcomp = VORtrainSeq('tTrain',t_inc, 'fps',[fpNS fpNorm]);
withpreobj = VORtrainSeqDiff('VORrel',VORrel, 'VORcomp',VORcomp);

exptobj=VORexperiment('WT',WTobj,'nopre',nopreobj,'withpre',withpreobj,...
    'pooled',pooled,'withprestyle','-','noprelabel','Sham stim','withprelabel','CF stim');

end

