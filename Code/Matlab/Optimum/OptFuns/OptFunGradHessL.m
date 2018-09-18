function [ f,gr,hess ] = OptFunGradHessL( x,s,fp,w,varargin )
%OPTFUNL function and gradient for fmincon on Laplace transform
%   x = parameters (off diagonal matrix elements)
%   f = function to be minimised (-SNR(t))
%   G = df/dx
%   hess = d^2f/dx^2

[Wp,Wm]=Params2Mats(x);

ev=ones(1,length(Wp));
q=Wp-Wm;
W=Wm+fp*q;

Zinv=ones(length(Wp))-W;
Zinvs=s*eye(length(Wp))+Zinv;

p=ev/Zinv;

Zw=Zinvs\w;
a=q*Zw;

f=-2*fp*(1-fp)*p*a;

if nargout>=2

    Za=Zinv\a;
    c=(p*q)/Zinvs;


    %dA(s)/dq_(i,j)
    dAdq = ((Zinvs\w)*p)';
    dAdq=dAdq-diag(dAdq)*ev;
    %dA(s)/dW^F_(i,j)
    dAdW = (Za*p)' + (Zw*c)';
    dAdW=dAdW-diag(dAdW)*ev;
    %dA/dWp_(i,i+1)+dA/dWm_(n-i+1,n-i)
    dAdWp=dAdq+fp*dAdW;
    dAdWm=-dAdq+(1-fp)*dAdW;

    gr = -Mats2Params(dAdWp,dAdWm);

end

if nargout>=3

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
    
    hess = Mats2ParamsHess(hesspp,hesspm,hessmp,hessmm);
    
end

    function tens=outer3(vec1,mat,vec2)
    %outer product of vector, matrix and vector
        tens=outer(vec1, outer(mat,vec2,false), true); 
    end

    function tens=subtractdiag(tens)
    %subtract diagonal element from each row for first two indices
        tens = tens - permute( outer( ev, sum(tens.*eyemask,1), true), [2 1 3 4]);
    end

    function tenstr=transposetens(tens)
    %swap first two and last two indices
        tenstr=permute(tens,[3 4 1 2]);
    end

    function [tens1,tens2]=subtracttranspose(tens)
    %subtract diagonal element from each row for first two indices, and
    %last two indices, and compute transposed tensor
        tens1 = subtractdiag(tens);
        tens2 = transposetens(tens1);
        tens2=subtractdiag(tens2);
        tens1=transposetens(tens2);
    end



% f=-real(crq(Wp,Wm,fp,w));
% f=-SNRarea( Wp, Wm, fp, w );
end

