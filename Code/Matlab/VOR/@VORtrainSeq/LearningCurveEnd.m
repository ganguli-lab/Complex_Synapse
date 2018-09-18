function [ S,Pt,t ] = LearningCurveEnd( obj,modelobj,dt )
%[S,p,t]=VORtrainSeq.LEARNINGCURVEEND(SynapseMemoryModel,dt) mean synaptic
%weight as function of time, starting in equilibrium state for
%VORtrainSeq.fps(1), then evolving according to VORtrainSeq.fps(2)...
%Only including final training epoch
%   dt = spacing of t values
%   S(t)=p(t)w
%   dp/dt = pW(fps(2))
%   p(0)W(fps(1))=0
%   actually uses max(t<t_n) in place of t_n.

error(CheckType(modelobj,'SynapseMemoryModel'));
error(CheckSize(modelobj,@isvalid));

tchanges = [0 obj.tTrain];
t = 0:dt:tchanges(end);

modelobj = modelobj.setFp(obj.fps(1));
%p0: ind(1,which state).
p0 = modelobj.EqProb;
%pt: ind(what time,which state).
pt = ones(length(t),1)*p0;
%S: ind(1,what time).
S = (pt*modelobj.w)';
Pt = pt;

for i=1:obj.numTrain
    %
    modelobj = obj.rs(i+1) * modelobj.setFp(obj.fps(i+1));

    valid = t>tchanges(i);
    ix = find(~valid,1,'last');
    
    p0 = pt(ix,:);%ind(1,which state).
    %V: ind(which state,which eigenmode).
    %D: ind(which eigenmode,which eigenmode)
    [V,D] = eig(modelobj.GetWf);
    expqt = exp(diag(D)*(t-t(ix)))';%ind(what time,which eigenmode).
    pt = (expqt*diag(p0*V))/V;%ind(what time,which state).
    
    newS = (pt*modelobj.w)';
    S(valid) = newS(valid);
    Pt(ix+1:end,:) = pt(ix+1:end,:);
end

valid = t>tchanges(end-1);
t(~valid) = [];
S(~valid) = [];
Pt(~valid,:) = [];

t = t - t(1);

end

