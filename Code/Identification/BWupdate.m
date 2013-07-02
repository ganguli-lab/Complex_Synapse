function [ M_new,initial_new,pstate ] = BWupdate( readouts,initial,outProj,M,potdep )
%[M_new,initial_new,pstate]=BWUPDATE(readouts,initial,outProj,M,potdep)
%Baum-Welch update of estiamted HMM
%   M_new       = updated Markov matrix/cell {Mpot,Mdep}
%   initial_new = updated prob dist of iniitial state (row vec)
%   pstate      = posterior prob of HMM being in each state at each time
%   readouts = which output was seen before each time-step 
%   initial  = prob dist of iniitial state (row vec)
%   outProj  = cell of diagonal matrices for each possible value of
%               output, with elements equal to prob of output
%either
%   M        = cell {Mpot,Mdep} of Markov matrices
%   potdep   = whether each transition was potentiating(1)/depressing(0)
%or
%   M        = Markov matrix

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

pstate=zeros(length(initial),length(readouts));
alpha=BWalpha(length(readouts),readouts,initial,outProj,M,potdep);
beta=BWbeta(1,readouts,outProj,M,potdep);
M_new={zeros(length(M{1}))};
if numel(M)==2
    M_new{2}=M_new{1};
end

for t=1:length(readouts)-1
    pstate(:,t)=alpha(t,:)'.*beta(:,t);
    pstate(:,t)=pstate(:,t)/sum(pstate(:,t));
    M_new{2-potdep(t)}=M_new{2-potdep(t)} + (beta(:,t+1)*alpha(t,:))' .* (M{2-potdep(t)}*outProj{readouts(t+1)});
end
    pstate(:,end)=alpha(end,:)'.*beta(:,end);
    pstate(:,end)=pstate(:,t)./sum(pstate(:,t));

initial_new=pstate(:,1)';
M_new{1}=StochastifyD(M_new{1});
if numel(M)==2
    M_new{2}=StochastifyD(M_new{2});
end



end

