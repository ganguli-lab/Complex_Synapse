function [ states,potdep ] = StateSequence( initial,M,randno,fp )
%[states,potdep]=STATESEQUENCE(initial,M,randno,fp) Monte-Carlo simulation
%of synaptic Markov process
%   states = sequance of occupied states
%   potdep = whether each transition was potentiating(1)/depressing(0)
%   initial = prob dist of iniitial state (row vec)
%either
%   M       = cell {Mpot,Mdep} of Markov matrices
%   randno  = matrix(2,n) of random numbers in [0,1]
%               first row controls whether transition is pot/dep
%               second row controls which transition is used
%   fp      =  fraction of transitions that are potentiating
%or
%   M       = Markov matrix
%   randno  = matrix(1,n) of random numbers in [0,1], controls which
%               transition is used .

error(CheckSize(initial,@isprob));
error(CheckValue(randno,@(x) all(all(inrange(x,0,1))),'inrange(0,1)'));
if iscell(M)
    error(CheckSize(M,@(x) numel(x)==2,'2 elements'));
    Mpot=M{1};
    Mdep=M{2};
    error(CheckSize(Mpot,@isstochasticD));
    error(CheckSize(Mdep,@isstochasticD));
    error(CheckSize(Mdep,@(x) samesize(x,Mpot),'samesize(Mpot)'));
    error(CheckSize(randno,@(x) size(x,1)==2,'2 rows'));
    error(CheckSize(fp,@isscalar));
    error(CheckValue(fp,@(x) inrange(x,0,1),'inrange(0,1)'));
    %
    Mpot=cumsum(Mpot,2);%change pdf -> cdf
    Mdep=cumsum(Mdep,2);%change pdf -> cdf
    potdep=(randno(1,1:end-1) < fp);
    randno=randno(2,:);
else
    error(CheckValue(M,@isstochasticD));
    error(CheckSize(randno,@isrow));
    %
    Mpot=cumsum(M,2);
    potdep=ones(1,size(randno,2)-1);
end
error(CheckSize(initial,@(x) length(x)==length(Mpot),'samesize(Mpot)'));
initial=cumsum(initial);%change pdf -> cdf

states=zeros(1,size(randno,2));

% states(1)=WhichBin(initial,randno(1));
states(1)=find(initial>randno(1),1,'first');

for i=2:length(states)
    if potdep(i-1)
        states(i)=find(Mpot(states(i-1),:)>randno(i),1,'first');
    else
        states(i)=find(Mdep(states(i-1),:)>randno(i),1,'first');
    end
end


end

