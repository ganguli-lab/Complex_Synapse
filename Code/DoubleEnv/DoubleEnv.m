function [ S1,whichcase,numexp ] = DoubleEnv( t1,t2,S2,n,varargin )
%DOUBLEENV (t1,t2,S2,n,...) Two-time envelope
%   t1 = time of maximisation
%   t2 = time of fixed SNR
%   S2 = SNR(t2)
%   n = # states
%
%   Case# = sum_i Theta(mu_i) 2^i-1
%   Constraints:
%       1: sum_a c_a q_a < 1
%       2: sum_a c_a     < n-1
%       3: c_a sqrt(q_a) < gamma
%   Parameter/Value pairs
%       t = time to evaluate SNR envelope (default=t1)
%       Constraint3 = Use constraint 3? (default=false)

t=t1;
Constraint3=false;
varargin=assignApplicable(varargin);

error(CheckSize(t,@isscalar));

if t1==t2
    numexp=1;
    whichcase=0;
    S1=S2*exp((t2-t)/t2);
    return;
end

if Constraint3
    numcases=8;
else
    numcases=4;
end%if Constraint3

if isscalar(t)
    siz=[numcases 2 1];
elseif isvector(t)
    siz=[numcases 2 length(t)];
else
    siz=[numcases 2 size(t)];
end
S1=zeros(siz);

for i=1:numcases
    for j=1:2
        S1(i,j,:) = feval(['Case' int2str(i-1) 'Exp' int2str(j)],...
            t1,t2,S2,n,'t',t,'Constraint3',Constraint3,varargin{:});
    end%for j
end%for i


[S1,numexp]=max(S1,[],2);
[S1,whichcase]=max(S1,[],1);
numexp=numexp(whichcase);
% S1=squeeze(S1);

end

