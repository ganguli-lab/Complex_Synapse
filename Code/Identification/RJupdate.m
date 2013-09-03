function [ M_new,initial_new,pstate,loglike ] = RJupdate( chunks,readouts,initial,outProj,M,potdep,varargin )
%[M_new,initial_new,pstate,loglike]=RJUPDATE(chunks,readouts,initial,outProj,M,potdep)
%Rabiner-Juang update of estiamted HMM
%   M_new       = updated Markov matrix/cell {Mpot,Mdep}
%   initial_new = updated prob dist of iniitial state (row vec)
%   pstate      = posterior prob of HMM being in each state at each time
%   loglike     = log likelihood of readouts under old model (prod over chunks)
%   chunks   = 2-by-K matrix of starts and ends of each chunk.
%   readouts = which output was seen before each time-step 
%   initial  = prob dist of iniitial state (row vec)
%   outProj  = cell of diagonal matrices for each possible value of
%               output, with elements equal to prob of output
%either
%   M        = cell {Mpot,Mdep} of Markov matrices
%   potdep   = whether each transition was potentiating(1)/depressing(0)
%or
%   M        = Markov matrix

lloffset=0;%avoid underflow by making this more negative
varargin=assignApplicable(varargin);

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
error(CheckSize(chunks,@ismatrix));
error(CheckSize(chunks,@(x) size(x,1)==2,'2-by-K'));
error(CheckValue(chunks,@(x) all(all(isint(x))),'2-by-K'));

pstate=zeros(length(initial),length(readouts));
M_new={zeros(length(M{1}))};
if numel(M)==2
    M_new{2}=M_new{1};
end
initial_new=zeros(size(initial));
loglike=0;

for i=1:size(chunks,2)
    range=chunks(1,i):chunks(2,i);
%    Weight=1/HMMlike(readouts(range),initial,outProj,M,potdep(range(1:end-1)));
    [chunkM,chunkInitial,chunkPs,chunkll]=BWupdate(readouts(range),initial,outProj,M,potdep(range(1:end-1)),'Normalise',false);
    Weight=exp(lloffset-chunkll);
    for j=1:length(chunkM)
        M_new{j}=M_new{j}+Weight*chunkM{j};
    end
    initial_new=initial_new+Weight*chunkInitial;
    pstate(:,range)=chunkPs*diag(1./sum(chunkPs,1));
    loglike=loglike+chunkll;
end
    
initial_new=initial_new/sum(initial_new);
M_new{1}=StochastifyD(M_new{1});
if numel(M)==2
    M_new{2}=StochastifyD(M_new{2});
else
    M_new=M_new{1};
end


end

