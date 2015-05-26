function [ newqv ] = PerturbC( qv,mode )
%PERTURBC Summary of this function goes here
%   Detailed explanation goes here

n=length(qv)/4;

switch mode
    case 's'
        %activity dependent plasticity
        qp=qv(1:n);
        qm=qv(n+1:2*n);
        [Wp,Wm]=MakeMultistate(qp,qm);
        %acticity independent plasticity
        qhp=qv(2*n+1:3*n);
        qhm=qv(3*n+1:end);
        [Hp,Hm]=MakeMultistate(qhp,qhm);
        H=Hp+Hm;

        p=(qp+qhp)./(qm+qhm);
        p=[1 cumprod(p)];
        p=p/sum(p);

        q=Wp-Wm;

        Zinv=ones(length(Wp))-Wm-0.5*q-H;
        Zinvs=s*eye(length(Wp))+Zinv;

        c=(p*q)/Zinvs;

        peturb=diff(c)<0;
    case 'c'
        q_pot_inc=qv(1:n);
        q_pot_dec=qv(n+1:2*n);
        q_dep_inc=qv(2*n+1:3*n);
        q_dep_dec=qv(3*n+1:end);

        Wp=StochastifyC(diag(q_pot_inc,1)+diag(q_pot_dec,-1));
        Wm=StochastifyC(diag(q_dep_inc,1)+diag(q_dep_dec,-1));

        p=(q_pot_inc+q_dep_inc)./(q_pot_dec+q_dep_dec);
        p=[1 cumprod(p)];
        p=p/sum(p);

        q=Wp-Wm;

        Zinv=ones(length(Wp))-Wm-0.5*q;
        Zinvs=s*eye(length(Wp))+Zinv;

         c=(p*q)/Zinvs;

        peturb=diff(c)<0;
end




end

