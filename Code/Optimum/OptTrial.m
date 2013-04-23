function [ lmax,KTp,KTm ] = OptTrial( nmax,n,varargin )
%[lmax,KTp,KTm]=OPTTRIAL(nmax,n) generate graphs for
%   lmax = is it a local max? Logical. Are all KTp,KTm >= 0? 
%   KTp,KTm = Kuhn-Tucker multipliers for maximisation of SNR(t) wrt elements of Wp,Wm
%   nmax = max #states
%   n = #states of trial uniform multistate model (default: nmax)

existsAndDefault('n',nmax);
Parent=gca;
UseDerivs=false;
InitRand=true;
varargin=assignApplicable(varargin);
varargin=[ {'Parent',Parent}, varargin];

cla(Parent);

t=10.^(-1:0.1:ceil(log10(nmax^2)));

[ h,yl ] = PlotEnvs( t,nmax,varargin{:} );
delete(h(3:end));

q=ones(1,n-1);
[Wp,Wm,w]=MakeSMS(q);

tm=TouchEnv(q);

s=SNRcurve(t,Wp,Wm,0.5,w);
plot(t,real(s),varargin{:});
smu=interp1(t,s,tm);

line([1;1]*tm,yl','Color','b','LineStyle','--',varargin{:});

[KTp,KTm]=KTmults(tm,Wp,Wm,0.5,w);

lmax=all(all( KTp>=0 & KTm>=0 ));

n=n+2;
tf=true;
while n>2 && tf
    n=n-2;
    w=[-ones(n/2,1);ones(n/2,1)];
    [Wp,Wm]=FindOpt(tm,n,'UseDerivs',UseDerivs,'InitRand',InitRand);
    tf=false;
%     [tf]=istransient(0.5*(Wp+Wm),1e-5,'UseP',true);
end
% if tf
%     Wp(ix,:)=[];
%     Wp(:,ix)=[];
%     Wm(ix,:)=[];
%     Wm(:,ix)=[];
%     w(ix)=[];
% end

s=SNRcurve(t,Wp,Wm,0.5,w);
plot(t,real(s),'r',varargin{:});
smn=interp1(t,s,tm);

% if smn>smu
    assignin('base','Wp',Wp);
    assignin('base','Wm',Wm);
%     disp(Wp);
%     disp(Wm);
% end



end

