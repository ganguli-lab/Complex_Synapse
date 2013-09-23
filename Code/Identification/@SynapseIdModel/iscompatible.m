function [ tf ] = iscompatible( modelobj,simobj )
%ISCOMPATIBLE are SynapseIdModel and SynapsePlastSeq compatible?
%   are indices in correct range?

tf = min(simobj.potdep)>0 &&...
    min(simobj.stateseq)>0 &&...
    min(simobj.readouts)>0 &&...
    max(simobj.potdep)<=length(modelobj.M) &&...
    max(simobj.stateseq)<=length(modelobj.initial) &&...
    max(simobj.readouts)<=length(modelobj.outProj) ;
    


end

