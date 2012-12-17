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

error(CheckSize(t1,@isscalar))
error(CheckValue(t1,@(x)x>0))
error(CheckSize(t2,@isscalar))
error(CheckValue(t2,@(x)x>0))
error(CheckSize(S2,@isscalar))
error(CheckValue(S2,@(x)x>0))
error(CheckSize(n,@isscalar))
error(CheckValue(n,@(x)x>0))

t=[];
Constraint3=false;
varargin=assignApplicable(varargin);

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

S1=zeros(numcases,2);

for i=1:numcases
    for j=1:2
        S1(i,j) = feval(['DoubleEnv_Case' int2str(i-1) 'Exp' int2str(j)],...
            t1,t2,S2,n,'Constraint3',Constraint3,varargin{:});
    end%for j
end%for i


[S1,numexp]=max(S1,[],2);
[S1,whichcase]=max(S1,[],1);
numexp=numexp(whichcase);
whichcase=whichcase-1;
% S1=squeeze(S1);

if ~isempty(t)
    S1 = feval(['Case' int2str(whichcase) 'Exp' int2str(numexp)],...
        t1,t2,S2,n,'t',t,'Constraint3',Constraint3,varargin{:});
end

end

