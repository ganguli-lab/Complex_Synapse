function [ S1 ] = Case2Exp1( t1,t2,S2,n,varargin )
%S1=CASE2EXP1(t1,t2,S2,n,...) Two-time envelope in case 2 with 1
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

error(CheckSize(t1,@isscalar))
error(CheckValue(t1,@(x)x>0))
error(CheckSize(t2,@isscalar))
error(CheckValue(t2,@(x)x>0))
error(CheckSize(S2,@isscalar))
error(CheckValue(S2,@(x)x>0))
error(CheckSize(n,@isscalar))
error(CheckValue(n,@(x)x>0))

gammasq=128/pi^4;

t=t1;
Constraint3=false;
varargin=assignApplicable(varargin);

c=n-1;
q=-lambertw(-S2*t2/(n-1))/t2;

S1=c*q*exp(-q*t);

valid = c*q<=1;
valid = valid && (t1-t2)/(1-q*t2)>=0;
if Constraint3
    valid = valid && c^2*q <= gammasq;
end%if

S1=S1*valid;

end%function

