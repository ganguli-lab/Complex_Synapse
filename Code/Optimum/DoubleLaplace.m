function [ A,dAdWp,dAdWm,p,fail ] = DoubleLaplace( s,Wp,Wm,fp,w )
%DOUBLELAPLACE Summary of this function goes here
%   Detailed explanation goes here




ev=ones(1,length(Wp));
q=Wp-Wm;
W=Wm+fp*q;

[U,S,~]=svd(W);

[smin,ix]=min(diag(S));
U=U(:,ix);
ix=U>0;

fail=smin < 1e-10 && (all(ix) || all(~ix));

if smin < 1e-10 && ~all(ix) && ~all(~ix)
    

    [A1,dAdWp1,dAdWm1,p1]=DoubleLaplace(s,Wp(ix,ix),Wm(ix,ix),fp,w(ix));
    [A2,dAdWp2,dAdWm2,p2]=DoubleLaplace(s,Wp(~ix,~ix),Wm(~ix,~ix),fp,w(~ix));

    f12=sum(p1*W(ix,~ix));
    f21=sum(p2*W(~ix,ix));

    f1=f21/(f12+f21);
    f2=f12/(f12+f21);

    dAdWp=zeros(size(W));
    dAdWm=zeros(size(W));
    p=zeros(1,length(W));

    A = (f1*A1 + f2*A2);

    p(ix) = f1*p1;
    p(~ix) = f2*p2;

    dAdWp(ix,ix) = f1*dAdWp1;
    dAdWp(~ix,~ix) = f2*dAdWp2;
    dAdWm(ix,ix) = f1*dAdWm1;
    dAdWm(~ix,~ix) = f2*dAdWm2;



    dAdWp(ix,~ix) = fp * (p1'*ones(1,sum(~ix))) * (-f1*A1 + f2*A2)/(f12+f21);
    dAdWp(~ix,ix) = fp * (p2'*ones(1,sum(ix))) * (f1*A1 - f2*A2)/(f12+f21);

    dAdWm(ix,~ix) = (1-fp) * (p1'*ones(1,sum(~ix))) * (-f1*A1 + f2*A2)/(f12+f21);
    dAdWm(~ix,ix) = (1-fp) * (p2'*ones(1,sum(ix))) * (f1*A1 - f2*A2)/(f12+f21);

else

    p=EqProb(W);

    Zinv=ones(length(Wp))-W;
    Zinvs=s*eye(length(Wp))+Zinv;

    a=q*(Zinvs\w);
    c=(p*q)/Zinvs;

    A=2*fp*(1-fp)*p*a;

    %dA(s)/dq_(i,j)
    dAdq = ((Zinvs\w)*p)';
    dAdq=dAdq-diag(dAdq)*ev;
    %dA(s)/dW^F_(i,j)
    dAdW = ((Zinv\a)*p)' + ((Zinvs\w)*c)';
    dAdW=dAdW-diag(dAdW)*ev;
    %dA/dWp_(i,i+1)+dA/dWm_(n-i+1,n-i)
    dAdWp=dAdq+fp*dAdW;
    dAdWm=-dAdq+(1-fp)*dAdW;

end








end

