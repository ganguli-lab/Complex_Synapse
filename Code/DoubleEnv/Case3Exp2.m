function [ S1 ] = Case3Exp2( t1,t2,S2,n,varargin )
%S1=CASE3EXP2(t1,t2,S2,n,...) Two-time envelope in case 1 with 2
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
varargin=assignApplicable(varargin);


q1guess=(1+exp(1)*t1/t2)/t2;
q2guess=1/t1;

qc=fsolve(@eqs,[q1guess,q2guess]);

q1=qc(1);
q2=qc(2);
c1=(1-(n-1)*q2)/(q1-q2);
c2=(1-(n-1)*q1)/(q2-q1);

S1=c1*q1*exp(-q1*t)+c2*q2*exp(-q2*t);

valid = mu1(q1,lambda(q1,q2))>=0;
valid = valid && mu2(q1,lambda(q1,q2))>=0;
if Constraint3
    valid = valid && c1^2*q1 <= gammasq;
    valid = valid && c2^2*q2 <= gammasq;
end%if

S1=S1*valid;



    function lam=lambda(qa,qb)
        lam=-(qa^2*t1*exp(-qa*t1)-qb^2*t1*exp(-qb*t1))...
            /(qa^2*t2*exp(-qa*t2)-qb^2*t2*exp(-qb*t2));
    end


    function mu=mu1(q,lmb)
        mu=(1-q*t1)*exp(-q*t1)+lmb*(1-q*t2)*exp(-q*t2);
    end

    function mu=mu2(q,lmb)
        mu=q^2*t1*exp(-q*t1)+lmb*q^2*t2*exp(-q*t2);
    end

    function val=eq1(q1,q2)
%         val=mu1(q1,lambda(q1,q2))-mu1(q2,lambda(q1,q2));
        val=(1-q1*t1)*t2*exp(-q1*t1) + (1-q1*t2)*t1*exp(-q1*t2)...
            - (1-q2*t1)*t2*exp(-q2*t1) - (1-q2*t2)*t1*exp(-q2*t2);
    end

%     function val=eq2(q1,q2,lambda)
%         val=mu2(q1,lambda)-mu2(q2,lambda);
%     end

    function val=eq2(q1,q2)
        val=(1-(n-1)*q2)*q1*exp(-q1*t2)-(1-(n-1)*q1)*q2*exp(-q2*t2)...
            - S2*(q1-q2);
    end

    function vals=eqs(qcs)
        vals=[0 0];
        vals(1)=eq1(qcs(1),qcs(2));
        vals(2)=eq2(qcs(1),qcs(2));
    end

end%function

