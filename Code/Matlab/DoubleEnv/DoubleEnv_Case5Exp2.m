function [ S1 ] = DoubleEnv_Case5Exp2( t1,t2,S2,n,varargin )
%S1=CASE7EXP2(t1,t2,S2,n,...) Two-time envelope in case 3 with 2
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

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='DoubleEnv_Case5Exp2';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addRequired('t1',@(x)validateattributes(x,{'numeric'},{'scalar','positive'},'DoubleEnv_Case5Exp2','t1',1));
    p.addRequired('t2',@(x)validateattributes(x,{'numeric'},{'scalar','positive'},'DoubleEnv_Case5Exp2','t2',2));
    p.addRequired('S2',@(x)validateattributes(x,{'numeric'},{'scalar','positive'},'DoubleEnv_Case5Exp2','S2',3));
    p.addRequired('n',@(x)validateattributes(x,{'numeric'},{'scalar','positive'},'DoubleEnv_Case5Exp2','n',4));
    p.addParameter('t',[],@(x) validateattributes(x,{'numeric'},{},'DoubleEnv_Case5Exp2','t'));
    p.addParameter('Constraint3',true,@(x) validateattributes(x,{'logical'},{'scalar'},'DoubleEnv_Case5Exp2','Constraint3'));
    p.addParameter('Display','off',@(x) parsevalidatestring(x,{'on','off'},'DoubleEnv_Case5Exp2','Display'));
    p.addParameter('Jacobian','on',@(x) parsevalidatestring(x,{'on','off'},'DoubleEnv_Case5Exp2','Jacobian'));
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



q1guess=0.9/t1;
q2guess=1.1/t2;

[qc,~,ef]=fsolve(@eqs,[q1guess,q2guess],optimset('Display',r.Display,'Jacobian',r.Jacobian));

q1=qc(1);
q2=qc(2);
c1=(1-(n-1)*q2)/(q1-q2);
c2=(1-(n-1)*q1)/(q2-q1);

S1=c1*q1*exp(-q1*r.t)+c2*q2*exp(-q2*r.t);

valid = r.Constraint3 && c1*q1+c2*q2<=1;
valid = valid && c1+c1<=n-1;
if r.Constraint3
    valid = valid && mu3(q1)>=0;
    valid = valid && mu3(q2) <= gammasq;
end%if

S1=S1*valid;

if ef ~=1
    S1=nan;
end




    function mu=mu3(q)
        mu=(t1-t2)/(1-2*q*t2);
    end

    function [val,grad]=eq1(q1,q2)
%         val=mu1(q1,lambda(q1,q2))-mu1(q2,lambda(q1,q2));
        val=0;
        grad(1)=-t1*(3-2*q1*t1)*(1-2*q2*t2)*exp(-q1*t1-q2*t2)...
            +-t2*(3-2*q1*t2)*(1-2*q2*t1)*exp(-q1*t2-q2*t1);
        grad(2)=-t2*(1-2*q1*t1)*(3-2*q2*t2)*exp(-q1*t1-q2*t2)...
            +t1*(1-2*q1*t2)*(3-2*q2*t1)*exp(-q1*t2-q2*t1);
    end

%     function val=eq2(q1,q2,lambda)
%         val=mu2(q1,lambda)-mu2(q2,lambda);
%     end

    function [val,grad]=eq2(q1,q2)
        val=0;
        grad(1)=sqrt(gammasq/q1)/2*(1-2*q1*t2)*exp(-q1*t2);
        grad(2)=sqrt(gammasq/q2)/2*(1-2*q2*t2)*exp(-q2*t2);
    end

    function [vals,jac]=eqs(qcs)
        vals=[0 0];
        [vals(1),grad1]=eq1(qcs(1),qcs(2));
        [vals(2),grad2]=eq2(qcs(1),qcs(2));
        jac=[grad1;grad2];
        vals(1)=vals(1)/(qcs(1)-qcs(2));
        jac(1,:)=(jac(1,:) + vals(1)*[-1,1])/(qcs(1)-qcs(2));
    end

end%function

