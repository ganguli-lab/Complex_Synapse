function [ lmax,KTp,KTm ] = OptTrial( nmax,varargin )
%[lmax,KTp,KTm]=OPTTRIAL(nmax,n) generate graphs for
%   lmax = is it a local max? Logical. Are all KTp,KTm >= 0? 
%   KTp,KTm = Kuhn-Tucker multipliers for maximisation of SNR(t) wrt elements of Wp,Wm
%   nmax = max #states
%   n = #states of trial uniform multistate model (default: nmax)

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='FindOpt';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addOptional('n',1,@(x)validateattributes(x,{'numeric'},{'scalar'},'FindOpt','reps',3));
    p.addParameter('Parent',gca,@(x) validateattributes(x,{'numeric'},{'scalar'},'FindOpt','InitRand'));
    p.addParameter('LineWidth',1,@(x) validateattributes(x,{'numeric'},{'scalar'},'FindOpt','InitRand'));
end
p.parse(varargin{:});
r=p.Results;
if any(strcmp('n',p.UsingDefaults))
    r.n=nmax;
end
plotArgs={'Parent',r.Parent,'LineWidth',r.LineWidth};
n=r.n;

cla(r.Parent);

t=10.^(-1:0.1:ceil(log10(nmax^2)));

[ h,yl ] = PlotEnvs( t,nmax,plotArgs{:} );
delete(h(3:end-1));

q=ones(1,n-1);
[Wp,Wm,w]=MakeSMS(q);

tm=TouchEnv(q);

s=SNRcurve(t,Wp,Wm,0.5,w);
plot(t,real(s),plotArgs{:});
% smu=interp1(t,s,tm);

line([1;1]*tm,yl','Color','b','LineStyle','--',plotArgs{:});

[KTp,KTm]=KTmultsSNR(tm,Wp,Wm,0.5,w);

lmax=all(all( KTp>=0 & KTm>=0 ));

n=n+2;
tf=true;
while n>2 && tf
    n=n-2;
    w=[-ones(n/2,1);ones(n/2,1)];
    [Wp,Wm]=FindOpt(tm,n,p.Unmatched);
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
plot(t,real(s),'r',plotArgs{:});
% smn=interp1(t,s,tm);

% if smn>smu
    assignin('base','Wp',Wp);
    assignin('base','Wm',Wm);
%     disp(Wp);
%     disp(Wm);
% end



end

