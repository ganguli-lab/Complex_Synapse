function [ alpha ] = BWalpha( t,readouts,initial,outProj,M,potdep )
%alpha=BWALPHA(t,readouts,initial,outProj,M,potdep) forward variables for Baum-Welch algorithm
%   alpha    = forward variables
%   t        = time-step up to which we want forward variables
%   readouts = which output was seen before each time-step 
%   initial  = prob dist of iniitial state (row vec)
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
error(CheckSize(initial,@isprob));
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
error(CheckSize(initial,@(x) length(x)==length(M{1}),'samesize(M)'));
for i=1:numel(outProj)
    error(CheckSize(outProj{i},@(x) samesize(x,M{1}),'samesize(M)'));
end

if any(potdep==0)
    potdep=2-potdep;
end

alpha=zeros(t,length(initial));
alpha(1,:)=initial*outProj{readouts(1)};

for i=2:t
    alpha(i,:)=alpha(i-1,:)*M{potdep(i-1)}*outProj{readouts(i)};
end

end

