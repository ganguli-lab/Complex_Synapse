function [ tf ] = iscompatible( simobj,modelobj )
%ISCOMPATIBLE are SynapseIdModel and SynapsePlastSeq compatible?
%   are indices in correct range?

tf = iscompatible@SynapsePlastSeq(simobj,modelobj) &&...
    simobj.NumStates<=modelobj.NumStates ;
    


end

