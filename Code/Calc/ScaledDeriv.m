function [ Delta ] = ScaledDeriv( Wp,Wm,fp,w,row,col )
%Delta=ScaledDeriv(WP,WM,FP,w,ROW,COL,PM) derivative of area wrt element
%   WP = potentiation transition rates
%   WM = depression transition rates
%   FP = Fraction of potentiation transitions
%   w = Weights of states (+/-1)
%   ROW,COL = which element we're changing
%   PM = +/-1=sgn(col-row), if we're changing WP/WM

pm=sign(col-row);
range=row:pm:col;

range=range(2:end-1);

W=fp*Wp+(1-fp)*Wm;
c=AreaCoeff(Wp,Wm,fp);
p=EqProb(W);
deta=DeltaEta(W,w);


deta=-diff(deta);
c=c-c(row);

p=p(row);
deta=deta(range);
c=c(range);

Delta=2*(0.5+pm*(fp-0.5))*p*c*deta;

end

