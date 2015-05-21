function [ zinv,zinvs ] = Zrconds( s,qv,mode )
%ZRCONDS Summary of this function goes here
%   Detailed explanation goes here


switch mode
    case 'n'
        qp=qv(1:length(qv)/2);
        qm=qv(length(qv)/2+1:end);
        [Wp,Wm]=MakeMultistate(qp,qm);

        q=Wp-Wm;

        Zinv=ones(length(Wp))-Wm-0.5*q;
        Zinvs=s*eye(length(Wp))+Zinv;
    case 's'
        n=length(qv)/4;

        %activity dependent plasticity
        qp=qv(1:n);
        qm=qv(n+1:2*n);
        [Wp,Wm]=MakeMultistate(qp,qm);
        %acticity independent plasticity
        qhp=qv(2*n+1:3*n);
        qhm=qv(3*n+1:end);
        [Hp,Hm]=MakeMultistate(qhp,qhm);
        H=Hp+Hm;

        q=Wp-Wm;

        Zinv=ones(length(Wp))-Wm-0.5*q-H;
        Zinvs=s*eye(length(Wp))+Zinv;
    case 'c'
        n=length(qv)/4;

        q_pot_inc=qv(1:n);
        q_pot_dec=qv(n+1:2*n);
        q_dep_inc=qv(2*n+1:3*n);
        q_dep_dec=qv(3*n+1:end);

        Wp=StochastifyC(diag(q_pot_inc,1)+diag(q_pot_dec,-1));
        Wm=StochastifyC(diag(q_dep_inc,1)+diag(q_dep_dec,-1));

        q=Wp-Wm;

        Zinv=ones(length(Wp))-Wm-0.5*q;
        Zinvs=s*eye(length(Wp))+Zinv;
    case 'a'
        n=length(qv)/2;

        q_pot_inc=qv(1:n);
        q_dep_dec=qv(n+1:2*n);
        q_pot_dec=(q_dep_dec-1).*(q_dep_dec>1);
        q_dep_inc=(q_pot_inc-1).*(q_pot_inc>1);

        Wp=StochastifyC(diag(q_pot_inc,1)+diag(q_pot_dec,-1));
        Wm=StochastifyC(diag(q_dep_inc,1)+diag(q_dep_dec,-1));

        q=Wp-Wm;

        Zinv=ones(length(Wp))-Wm-0.5*q;
        Zinvs=s*eye(length(Wp))+Zinv;
end


zinv=rcond(Zinv);
zinvs=rcond(Zinvs);

end

