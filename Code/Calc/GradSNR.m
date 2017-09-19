function [ S, dSdWp,dSdWm ] = GradSNR( t, Wp, Wm, fp, w, varargin  )
%[S,GS]=GRADSNR(T,WP,WM,FP,w) gradient of SNR curve wrt matrix elements
%   T = time value
%   WP = potentiation transition rates
%   WM = depression transition rates
%   FP = Fraction of potentiation transitions
%   w = Weights of states (+/-1)

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='GradSNR';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addRequired('t',@(x)validateattributes(x,{'numeric'},{'scalar'},'GradSNR','t',1));
    p.addRequired('Wp',@(x)validateattributes(x,{'numeric'},{'2d','square'},'GradSNR','Wp',2));
    p.addRequired('Wm',@(x)validateattributes(x,{'numeric'},{'2d','square'},'GradSNR','Wm',3));
    p.addRequired('fp',@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1},'GradSNR','fp',4));
    p.addRequired('w',@(x)validateattributes(x,{'numeric'},{'column'},'GradSNR','w',5));
    p.addParameter('NoDiag',true,@(x)validateattributes(x,{'logical'},{'scalar'},'GradSNR','NoDiag'));
    p.addParameter('thresh',1e-7,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'},'GradSNR','thresh'));
    p.addParameter('DegThresh',1e-3,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'},'GradSNR','DegThresh'));
    p.addParameter('RCondThresh',1e-5,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'},'GradSNR','RCondThresh'));
    p.addParameter('UseZ',true,@(x)validateattributes(x,{'logical'},{'scalar'},'GradSNR','UseZ'));
    p.addParameter('UseV',false,@(x)validateattributes(x,{'logical'},{'scalar'},'GradSNR','UseV'));
end
p.parse(t,Wp,Wm,fp,w,varargin{:});
r=p.Results;

error(CheckSize(r.Wm,@(x)samesize(r.Wp,x),'samesize(Wp)'));%also square matrix of same size
error(CheckSize(r.w,@(x)length(x)==length(r.Wp),'samesize(Wp)'));


q=r.Wp-r.Wm;
W=r.fp*q + r.Wm;

[u,qa]=eig(-W);
qa=diag(qa);
[qa,ix]=sort(qa);
u=u(:,ix);

Zinv=ones(length(r.w)) - W;

%this method doesn't work when eigenvalues are nearly degenerate
if min(diff(qa)) > r.DegThresh
    [v,qb]=eig(-W');
    qb=diag(qb);
    [~,ix]=sort(qb);
    v=v(:,ix);
    v=diag(1./diag(v.'*u)) * v.';

    Z=u * diag(1./[1;qa(2:end)]) * v;
    pr=v(1,:);
    pr=pr/sum(pr);
    r.UseV=true;
%this method doesn't work when eigenvectors are nearly parallel or W is
%nearly non-ergodic
elseif rcond(Zinv) > r.RCondThresh
    pr = ones(1,length(r.w))/Zinv;
    r.UseZ=false;
%this method doesn't work when eigenvectors are nearly parallel
else
    Z=u * diag(1./[1;qa(2:end)]) / u;
    pr=[1 zeros(1,length(qa)-1)] / u;
    pr=pr/sum(pr);
end

expLt=exp(-qa*t);
if r.UseV
    expWt=u*diag(expLt)*v;
else
    expWt=u*diag(expLt)/u;
end

S=pr*q*expWt*r.w;

%deriv wrt q_ij
dSdq=((expWt*r.w)*pr).';

%deriv wrt W_ij, due to p
if r.UseZ
    dSdp=((Z*q*expWt*r.w)*pr).';
else
    dSdp=((Zinv\q*expWt*r.w)*pr).';
end

%deriv wrt W_ij, due to expWt
%ref: http://dx.doi.org/10.1002/nme.263
FF=expLt*ones(1,length(r.w));
F=FF-FF.';
qa=qa*ones(1,length(r.w));
qa=qa.'-qa;
%check for degenerate evals
degenerate=abs(qa)<r.thresh;
qa(degenerate)=1;
F(degenerate)=t*FF(degenerate);
%
F=F./qa;
if r.UseV
    dSdexpWt=(u*diag(v*r.w)*F*diag(pr*q*u)*v).';
else
    dSdexpWt=(u*diag(u\r.w)*F*diag(pr*q*u)/u).';
end

%deriv wrt W_ij
dSdW=dSdp+dSdexpWt;

if r.NoDiag
    dSdq=dSdq-diag(dSdq)*ones(1,length(r.w));
    dSdW=dSdW-diag(dSdW)*ones(1,length(r.w));
end

dSdWp=r.fp*dSdW+dSdq;
dSdWm=dSdW-dSdWp;

dSdWp=real(dSdWp);
dSdWm=real(dSdWm);

end

