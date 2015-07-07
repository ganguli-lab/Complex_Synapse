function [hesspp,hesspm,hessmp,hessmm]=SNRLaplaceHess(s,Wp,Wm,fp,w)
%[hesspp,hesspm,hessmp,hessmm]=SNRLAPLACEHESS(s,Wp,Wm,fp,w) hessian of
%Laplace transform of SNR curve

eyemask=outer(eye(length(Wp)),ones(length(Wp)),false);

fm=1-fp;

ev=ones(1,length(Wp));
q=Wp-Wm;
W=Wm+fp*q;

Zinv=ones(length(Wp))-W;

p=ev/Zinv;

Zinvs=s*eye(length(Wp))+Zinv;

Zw=Zinvs\w;
a=q*Zw;
Za=Zinv\a;
c=(p*q)/Zinvs;

Z=inv(Zinv);
Zs=inv(Zinvs);

ZqZ = (Zinv \ q) / Zinvs;

hess2 = outer3(p,Z,Za); 
[hess2,hess1]=subtracttranspose(hess2);

hess4 = outer3(p,ZqZ,Zw);
[hess4,hess7]=subtracttranspose(hess4);

hess10 = outer3(c,Zs,Zw);
[hess10,hess9]=subtracttranspose(hess10);

hess3 = outer(p,Z,Zw);
[hess3,hess5]=subtracttranspose(hess3);

hess6 = outer(p,Zs,Zw);
[hess6,hess8]=subtracttranspose(hess6);

hessWW = hess1 + hess2 + hess4 + hess7 + hess9 + hess10;
hessWq = hess3 + hess8;
hessqW = hess5 + hess6;

hesspp = fp^2 * hessWW + fp * ( hessWq + hessqW );
hesspm = fp*fm * hessWW - fp * hessWq + fm * hessqW;
hessmp = fp*fm * hessWW + fm * hessWq - fp * hessqW;
hessmm = fm^2 * hessWW - fm * ( hessWq + hessqW );

    function tens=outer3(vec1,mat,vec2)
        tens=outer(vec1, outer(mat,vec2,true), true); 
    end

    function tens=subtractdiag(tens)
        tens = tens - outer( ev, sum(tens.*eyemask,1), true);
    end

    function tenstr=transposetens(tens)
        tenstr=permute(tens,[3 4 1 2]);
    end

    function [tens1,tens2]=subtracttranspose(tens)
        tens1 = subtractdiag(tens);
        tens2 = transposetens(tens1);
        tens2=subtractdiag(tens2);
        tens1=transposetens(tens2);
    end

end





