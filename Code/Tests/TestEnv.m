function TestEnv(trialtype,varargin)

existsAndDefault('trialtype','numulti');

randtrans=false;
pmulti=false;
numulti=false;

eval([trialtype '=true;']);
if ~(randtrans||pmulti||numulti)
    numulti=true;
end

fp=0.5;
len=20;
numtrials=100;
if randtrans || pmulti
    sp=0.31;
end%if
if pmulti || numulti
    di=0.7;
    sh=0.1;
end%if

varargin=assignApplicable(varargin);

if randtrans
    w=[-ones(len/2,1);ones(len/2,1)];
end%if
if numulti
    qu=di*ones(1,len-1);
%     numtrials=len-1;
end%if

senv=[];
numfail=0;
for i=1:numtrials
    %generate trial
    %random
    if randtrans
        W=RandTrans(len,sp);
        [Wp,Wm,neww]=TriangleDcmp(W,fp,w);
        Wpu=Wp;
        Wmu=Wm;
    end%if
    %perturb uniform multistate
    if pmulti
        [Wp,Wm,neww]=SMSneighbour(len,fp,di,sh,sp);
        [Wpu,Wmu]=SMSneighbour(len,fp,di,0,sp);
    %non-uniform multistate;
    end%if
    if numulti
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
    if rcond(fp*Wp+(1-fp)*Wm-ones(size(Wp))) < 1e-15
        continue;
    end
    [tf,s,t]=CheckEnv(Wp,Wm,fp,neww);
    [env,bnds]=SNRenvelope(t,len);
    [~,su]=CheckEnv(Wpu,Wmu,fp,oldw);
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
        bnds2=len*(len-2)/(4*di);
        bnds2=[bnds2;len*(len+2)/(4*di)];
        bnds3=len^2/pi^2;
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