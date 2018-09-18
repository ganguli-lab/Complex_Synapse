function TestEnv(varargin)

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='TestEnv';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addOptional('trialtype','numulti',@(x)parsevalidatestring(x,{'nmulti','pmulti','randtrans'},'TestEnv','thresh',1));
    p.addParameter('fp',0.5,@(x)validateattributes(x,{'numeric'},{'scalar'},'TestEnv','fp'));
    p.addParameter('len',20,@(x)validateattributes(x,{'numeric'},{'scalar','integer'},'TestEnv','len'));
    p.addParameter('numtrials',100,@(x)validateattributes(x,{'numeric'},{'scalar','integer'},'TestEnv','numtrials'));
    p.addParameter('sp',0.31,@(x)validateattributes(x,{'numeric'},{'scalar'},'TestEnv','sp'));
    p.addParameter('di',0.7,@(x)validateattributes(x,{'numeric'},{'scalar'},'TestEnv','di'));
    p.addParameter('sh',0.1,@(x)validateattributes(x,{'numeric'},{'scalar'},'TestEnv','sh'));
end
p.parse(varargin{:});
r=p.Results;

S.randtrans=false;
S.pmulti=false;
S.numulti=false;

S.(r.trialtype)=true;

% fp=0.5;
% len=20;
% numtrials=100;
% if randtrans || pmulti
%     sp=0.31;
% end%if
% if pmulti || numulti
%     di=0.7;
%     sh=0.1;
% end%if
% 

if S.randtrans
    w=[-ones(r.len/2,1);ones(r.len/2,1)];
end%if
if r.numulti
    qu=r.di*ones(1,r.len-1);
%     numtrials=len-1;
end%if

senv=[];
numfail=0;
for i=1:r.numtrials
    %generate trial
    %random
    if S.randtrans
        W=RandTrans(r.len,r.sp);
        [Wp,Wm,neww]=TriangleDcmp(W,r.fp,w);
        Wpu=Wp;
        Wmu=Wm;
    end%if
    %perturb uniform multistate
    if S.pmulti
        [Wp,Wm,neww]=SMSneighbour(r.len,r.fp,r.di,r.sh,r.sp);
        [Wpu,Wmu]=SMSneighbour(r.len,r.fp,r.di,0,r.sp);
    %non-uniform multistate;
    end%if
    if S.numulti
        qp=qu;
        qm=qu;
        qp=qp(1:end-2);
        qm=qm(1:end-2);
%         qp(len/2)=sh+(di-sh)*(numtrials-i)/(numtrials-1);
%         qm(len/2)=qp(len/2);
%         qp(len-1)=sh+(di-sh)*(numtrials-i)/(numtrials-1);
%         qm(1)=qp(len/2);
%         qp(i)=sh;
%         qm(i)=sh;
%         qp(len-i)=sh;
%         qm(len-i)=sh;
        [Wp,Wm,neww]=MakeMultistate(qp,qm);
        [Wpu,Wmu,oldw]=MakeMultistate(qu,qu);
    end%if
    %test trial
    if rcond(r.fp*Wp+(1-r.fp)*Wm-ones(size(Wp))) < 1e-15
        continue;
    end
    [tf,s,t]=CheckEnv(Wp,Wm,r.fp,neww);
    [env,bnds]=SNRenvelope(t,r.len);
    [~,su]=CheckEnv(Wpu,Wmu,r.fp,oldw);
    if isempty(senv)
        senv=real(s);
    else
        senv=max(senv,real(s));
    end
    if mod(i,10)==0 %|| numulti
        disp(int2str(i));
%         disp(num2str(qp(len/2)));
    end%if
    if tf || mod(i,10)==0 %|| numulti
        bnds2=r.len*(r.len-2)/(4*r.di);
        bnds2=[bnds2;r.len*(r.len+2)/(4*r.di)];
        bnds3=r.len^2/pi^2;
        yl=[env(end) 1];
        plot(t,abs(real([s;su(1:length(s));env(1:length(s))])))
        line([1;1]*bnds',yl'*ones(1,length(bnds)),'Color','k');
        line([1;1]*bnds2',yl'*ones(1,length(bnds2)),'Color','m');
        line([1;1]*bnds3',yl'*ones(1,length(bnds3)),'Color','c');
        set(gca,'XScale','log','YScale','log');
        ylim(yl);
        xlim([t(2) t(end)]);
        xlabel('t');
        ylabel('SNR');
        legend({'perturbed','unperturbed','envelope'},'Location','Best');
        pause;
    end%if
    if tf
        numfail=numfail+1;
    end%if
end%for i

yl=[env(end) 1];
plot(t,abs(real([senv;su;env])))
line([1;1]*bnds',yl'*ones(1,length(bnds)),'Color','k');
set(gca,'XScale','log','YScale','log');
ylim(yl);
xlim([t(2) t(end)]);
xlabel('t');
ylabel('SNR');
legend({'max perturbed','unperturbed','envelope'},'Location','Best');

disp(['number of failures = ' int2str(numfail)]);
end
% clear len sp fp di sh numfail Wp Wm neww tf i su W Wmu Wpu w bnds qp qm qpu qmu s t env;