function [ tau ] = MixTime( W )
%tau=MIXTIME(W) mixing time, inverse of -(max(Re(eigenvalue~=0)) 
%   Detailed explanation goes here

evs=sort(real(eig(-W)));
tau=1/evs(2);


end

