function [ Obs ] = MakeObs( modelobj )
%Obs=MAKEOBS(modelobj) construct observation matrix
%   Obs(i,x) = probability of observing x modelobj.when in state i

Obs=zeros(modelobj.NumStates,modelobj.NumWvals);

inds=sub2ind(size(Obs),1:modelobj.NumStates,modelobj.w');

Obs(inds)=1;


end

