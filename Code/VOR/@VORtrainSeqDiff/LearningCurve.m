function [ S,Pt,t,Pt_other ] = LearningCurve( obj,modelobj,dt )
%[S,p,t]=VORtrainSeqDiff.LEARNINGCURVE(SynapseMemoryModel,dt) mean synaptic
%weight as function of time, starting in equilibrium state for
%VORtrainSeq.fps(1), then evolving according to VORtrainSeq.fps(2)...
%   dt = spacing of t values
%   S(t)=p(t)w
%   dp/dt = pW(fps(2))
%   p(0)W(fps(1))=0
%   actually uses max(t<t_n) in place of t_n.

error(CheckType(modelobj,'SynapseMemoryModel'));
error(CheckSize(modelobj,@isvalid));

[S_rel,Pt,t] = obj.VORrel.LearningCurve(modelobj, dt);

[S_other, Pt_other] = obj.VORcomp.LearningCurve(modelobj, dt);

S = (1 - obj.frac_other) * S_rel - obj.frac_other * S_other;


end

