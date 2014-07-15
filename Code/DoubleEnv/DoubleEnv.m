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

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='DoubleEnv';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addRequired('t1',@(x)validateattributes(x,{'numeric'},{'scalar','positive'},'DoubleEnv','t1',1));
    p.addRequired('t2',@(x)validateattributes(x,{'numeric'},{'scalar','positive'},'DoubleEnv','t2',2));
    p.addRequired('S2',@(x)validateattributes(x,{'numeric'},{'scalar','positive'},'DoubleEnv','S2',3));
    p.addRequired('n',@(x)validateattributes(x,{'numeric'},{'scalar','positive'},'DoubleEnv','n',4));
    p.addParameter('t',[],@(x) validateattributes(x,{'numeric'},{},'DoubleEnv','t'));
    p.addParameter('Constraint3',false,@(x) validateattributes(x,{'logical'},{'scalar'},'DoubleEnv','Constraint3'));
end
p.parse(t1,t2,S2,n,varargin{:});
r=p.Results;

t=r.t;
t1=r.t1;
t2=r.t2;
S2=r.S2;

if t1==t2
    numexp=1;
    whichcase=0;
    if isempty(t)
        S1=S2;
    else
        S1=S2*exp((t2-t)/t2);
    end
    return;
end

if r.Constraint3
    numcases=8;
else
    numcases=4;
end%if Constraint3

S1=zeros(numcases,2);

for i=1:numcases
    for j=1:2
        S1(i,j) = feval(['DoubleEnv_Case' int2str(i-1) 'Exp' int2str(j)],...
            t1,t2,S2,r.n,'Constraint3',r.Constraint3,p.Unmatched);
    end%for j
end%for i


[S1,numexp]=max(S1,[],2);
[S1,whichcase]=max(S1,[],1);
numexp=numexp(whichcase);
whichcase=whichcase-1;
% S1=squeeze(S1);

if ~isempty(t)
    S1 = feval(['Case' int2str(whichcase) 'Exp' int2str(numexp)],...
        t1,t2,S2,r.n,'t',t,'Constraint3',r.Constraint3,p.Unmatched);
end

end

