function [ S1 ] = DoubleEnv_Case2Exp2( t1,t2,S2,n,varargin )
%S1=CASE2EXP2(t1,t2,S2,n,...) Two-time envelope in case 2 with 2
%exponentials
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
Display='off';
Jacobian='on';
varargin=assignApplicable(varargin);


q1guess=0.9/t1;
q2guess=1.1/t2;

[qc,~,ef]=fsolve(@eqs,[q1guess,q2guess],optimset('Display',Display,'Jacobian',Jacobian));

q1=qc(1);
q2=qc(2);
c1=(S2-(n-1)*q2*exp(-q2*t2))/(q1*exp(-q1*t2)-q2*exp(-q2*t2));
c2=n-1-c1;

S1=c1*q1*exp(-q1*t)+c2*q2*exp(-q2*t);

valid = c1*q1+c2*q2<=1;
valid = valid && mu2(q1)>=0;
if Constraint3
    valid = valid && c1^2*q1 <= gammasq;
    valid = valid && c2^2*q2 <= gammasq;
end%if

S1=S1*valid;

if ef ~=1
    S1=nan;
end

%     function lam=lambda(q)
%         lam=-(1-q*t1)/(1-q*t2)*exp(q*(t2-t1));
%     end

    function mu=mu2(q)
        mu=(t1-t2)/(1-q*t2);
    end

    function [val,grad]=eq1(q1,q2)
        val=(1-q1*t1)*(1-q2*t2)*exp(-q1*t1-q2*t2)-(1-q1*t2)*(1-q2*t1)*exp(-q1*t2-q2*t1);
        grad(1)=-t1*(2-q1*t1)*(1-q2*t2)*exp(-q1*t1-q2*t2)+t2*(2-q1*t2)*(1-q2*t1)*exp(-q1*t2-q2*t1);
        grad(2)=-t2*(1-q1*t1)*(2-q2*t2)*exp(-q1*t1-q2*t2)+t1*(1-q1*t2)*(2-q2*t1)*exp(-q1*t2-q2*t1);
    end

    function [val,grad]=eq2(q1,q2)
        val=(1-q2*t2)*q1^2*exp(-q1*t1)-(1-q1*t2)*q2^2*exp(-q2*t1);
        grad(1)=(1-q2*t2)*q1*(2-q1*t1)*exp(-q1*t1)+t2*q2^2*exp(-q2*t1);
        grad(2)=-t2*q1^2*exp(-q1*t1)-(1-q1*t2)*q2*(2-q2*t1)*exp(-q2*t1);
    end

    function [vals,jac]=eqs(qcs)
        vals=[0 0];
        [vals(1),grad1]=eq1(qcs(1),qcs(2));
        [vals(2),grad2]=eq2(qcs(1),qcs(2));
        vals=vals/(qcs(1)-qcs(2));
        jac=([grad1;grad2] + vals'*[-1,1])/(qcs(1)-qcs(2));
    end

end%function

