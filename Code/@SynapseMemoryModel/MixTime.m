function [ tau ] = MixTime( obj )
%tau=SynapseMemoryModel.MIXTIME() mixing time, inverse of -(max(Re(eigenvalue~=0)) 
%   Detailed explanation goes here

evs=sort(real(eig(-obj.GetWf)));
tau=1/evs(2);


end

