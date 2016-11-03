function [ exptobj ] = VORbuilderChR2( builder_h,numStates,paramPot,paramDep,fpNorm,fpInc,fpNS,t_inc,t_NS,pooled )
%VORexptobj=VORBUILDERCHR2(builder_h,numStates,paramPot,paramDep,fpNorm,fpInc,fpNS,t_inc,t_NS,pooled)
%building a VORexperiment object with non-specific ChR2 CF stim
%   builder_h  = function handle [Wp,Wm,w]=builder_h(numStates,param)
%   num_states = number of synaptic states
%   paramPot   = paramter for Wp
%   paramDep   = paramter for Wm
%   fpNorm     = fraction of potentiating events at baseline
%   fpInc      = fraction of potentiating events during gain-increase training
%   fpNS       = fraction of potentiating events during non-specific CF stim
%   t_inc      = duration of gain-increase training
%   t_NS      = duration of non-specific CF stim
%   pooled     = is it a pooled resource model?



[Wp,~,w] = builder_h(numStates,paramPot);
[~,Wm] = builder_h(numStates,paramDep);

WTobj = SynapseMemoryModel('Wp',Wp,'Wm',Wm,'w',w,'fp',fpNorm);

VORrel = VORtrainSeq('tTrain',[t_NS t_NS+t_inc], 'fps',[fpNorm fpNorm fpInc]);
VORcomp = VORtrainSeq('tTrain',[t_NS t_NS+t_inc], 'fps',[fpNorm fpNorm fpNorm]);
nopreobj = VORtrainSeqDiff('VORrel',VORrel, 'VORcomp',VORcomp);

VORrel = VORtrainSeq('tTrain',[t_NS t_NS+t_inc], 'fps',[fpNorm fpNS fpInc]);
VORcomp = VORtrainSeq('tTrain',[t_NS t_NS+t_inc], 'fps',[fpNorm fpNS fpNorm]);
withpreobj = VORtrainSeqDiff('VORrel',VORrel, 'VORcomp',VORcomp);

exptobj=VORexperiment('WT',WTobj,'nopre',nopreobj,'withpre',withpreobj,...
    'pooled',pooled,'withprestyle','-','noprelabel','Sham stim','withprelabel','CF stim');

end

