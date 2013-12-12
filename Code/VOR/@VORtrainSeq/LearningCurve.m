function [ S,Pt,t ] = LearningCurve( obj,modelobj )
%S=LEARNINGCURVE(W0,w,t,W1,...) mean synaptic weight as function of time,
%starting in equilibrium state for W0, then evolving according to W1
%   S(t)=p(t)w
%   dp/dt = pW
%   p(0)W_0=0
%   actually uses max(t<t_n) in place of t_n.

error(CheckType(modelobj,'SynapseMemoryModel'));
error(CheckSize(modelobj,@isvalid));

tchanges=[0 obj.tTrain];
t=0:tchanges(end);

modelobj=modelobj.SetFp(obj.fps(1));
%p0: ind(1,which state).
p0=modelobj.EqProb;
%pt: ind(what time,which state).
pt=ones(length(t),1)*p0;
%S: ind(1,what time).
S=(pt*modelobj.w)';
Pt=pt;

for i=1:obj.numTrain
    %
    modelobj=modelobj.SetFp(obj.fps(i+1));

    valid=t>tchanges(i);
    ix=find(~valid,1,'last');
    
    p0=pt(ix,:);%ind(1,which state).
    %V: ind(which state,which eigenmode).
    %D: ind(which eigenmode,which eigenmode)
    [V,D]=eig(modelobj.GetWf);
    expqt=exp(diag(D)*(t-t(ix)))';%ind(what time,which eigenmode).
    pt=(expqt*diag(p0*V))/V;%ind(what time,which state).
    
    newS=(pt*modelobj.w)';
    S(valid)=newS(valid);
    Pt(ix+1:end,:)=pt(ix+1:end,:);
end



end

