function [ beta ] = BWbeta( t,readouts,outProj,M,potdep )
%alpha=BWBETA(t,readouts,initial,outProj,M,potdep) forward variables for Baum-Welch algorithm
%   alpha    = forward variables
%   t        = time-step from which we want forward variables
%   readouts = which output was seen before each time-step 
%   outProj  = cell of diagonal matrices for each possible value of
%               output, with elements equal to prob of output
%either
%   M        = cell {Mpot,Mdep} of Markov matrices
%   potdep   = whether each transition was potentiating(1)/depressing(0)
%or
%   M        = Markov matrix


error(CheckSize(t,@isscalar));
error(CheckValue(t,@isint));
error(CheckSize(readouts,@isvector));
error(CheckValue(readouts,@(x) all(isint(x)),'isint'));
error(CheckType(outProj,'cell'));
error(CheckSize(outProj,@(x) numel(x)>=max(readouts),'numel>=max(readouts)'));
if iscell(M)
    error(CheckSize(M,@(x) numel(x)==2,'2 elements'));
    error(CheckSize(M{1},@isstochasticD));
    error(CheckSize(M{2},@isstochasticD));
    error(CheckSize(M{2},@(x) samesize(x,M{1}),'samesize(Mpot)'));
    error(CheckSize(potdep,@(x) samesize(x,readouts(2:end)),'samesize(readouts(2:end))'));
    %
else
    error(CheckValue(M,@isstochasticD));
    %
    M={M};
    potdep=ones(1,size(readouts,2)-1);
end
for i=1:numel(outProj)
    error(CheckSize(outProj{i},@(x) samesize(x,M{1}),'samesize(M)'));
end

beta=ones(length(M{1}),length(readouts)-t+1);

for i=length(readouts):-1:t+1
    beta(:,i-t)=M{2-potdep(i-1)}*outProj{readouts(i)}*beta(:,i-t+1);
end

end

