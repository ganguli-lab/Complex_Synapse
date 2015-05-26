function [ newqv ] = PerturbK( qv,mode )
%PERTURBC Summary of this function goes here
%   Detailed explanation goes here

n=length(qv)/4;

switch mode
    case 's'
        %activity dependent plasticity
        qp=qv(1:n);
        qm=qv(n+1:2*n);
        %acticity independent plasticity
        qhp=qv(2*n+1:3*n);
        qhm=qv(3*n+1:end);
        
        dqp=min(qhp,(1-qp)/2);
        dqm=min(qhm,(1-qm)/2);
        
        newqp=qp+2*dqp;
        newqhp=qhp-dqp;
        
        newqm=qm+2*dqm;
        
        newqv=[newqp newqm newqhp newqhm];
    case 'c'
        q_pot_inc=qv(1:n);
        q_pot_dec=qv(n+1:2*n);
        q_dep_inc=qv(2*n+1:3*n);
        q_dep_dec=qv(3*n+1:end);

end




end

