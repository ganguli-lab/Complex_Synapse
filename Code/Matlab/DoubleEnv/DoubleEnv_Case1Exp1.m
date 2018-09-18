function [ S1 ] = DoubleEnv_Case1Exp1( t1,t2,S2,n,varargin )
%S1=CASE1EXP1(t1,t2,S2,n,...) Two-time envelope in case 1 with 1
%exponential 
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
    p.FunctionName='DoubleEnv_Case1Exp1';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addRequired('t1',@(x)validateattributes(x,{'numeric'},{'scalar','positive'},'DoubleEnv_Case1Exp1','t1',1));
    p.addRequired('t2',@(x)validateattributes(x,{'numeric'},{'scalar','positive'},'DoubleEnv_Case1Exp1','t2',2));
    p.addRequired('S2',@(x)validateattributes(x,{'numeric'},{'scalar','positive'},'DoubleEnv_Case1Exp1','S2',3));
    p.addRequired('n',@(x)validateattributes(x,{'numeric'},{'scalar','positive'},'DoubleEnv_Case1Exp1','n',4));
    p.addParameter('t',[],@(x) validateattributes(x,{'numeric'},{},'DoubleEnv_Case1Exp1','t'));
    p.addParameter('Constraint3',false,@(x) validateattributes(x,{'logical'},{'scalar'},'DoubleEnv_Case1Exp1','Constraint3'));
end
p.parse(t1,t2,S2,n,varargin{:});
r=p.Results;
if any(strcmp('t',p.UsingDefaults))
    r.t=r.t1;
end

t1=r.t1;
t2=r.t2;
S2=r.S2;

gammasq=128/pi^4;



S1=S2.^(r.t/t2);

valid = t1<=t2;
valid = valid && t2 <= - (n-1)*log(S2);
if r.Constraint3
    valid = valid && t2 <= - gammasq*log(S2);
end%if

S1=S1*valid;

end%function

