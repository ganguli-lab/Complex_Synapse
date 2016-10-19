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

[S_rel,Pt,t] = LearningCurve@VORtrainSeq(obj, modelobj, dt);

fps_old = obj.fps;
obj = obj.setFp(obj.fps_other);

[S_other, Pt_other] = LearningCurve@VORtrainSeq(obj, modelobj, dt);

obj = obj.setFp(fps_old);

S = (1 - obj.frac_other) * S_rel - obj.frac_other * S_other;


end

