function [tmin,tmax,pert]=Crossover( varargin )
%CROSSOVER Summary of this function goes here
%   Detailed explanation goes here

fp=0.5;
len=20;
numtrials=100;
di=1;
sh=0.1;

varargin=assignApplicable(varargin);

qu=di*ones(1,len-1);

pert=zeros(1,numtrials);
tmin=pert;
tmax=pert;

for i=1:numtrials
    pert(i)=sh+(di-sh)*(numtrials-i)/(numtrials);
    %generate trial
        qp=qu;
        qm=qu;
        qp(len-1)=pert(i);
        qm(1)=pert(i);
        [Wp,Wm,neww]=MakeMultistate(qp,qm);
        [Wpu,Wmu]=MakeMultistate(qu,qu);
    %test trial
    if rcond(fp*Wp+(1-fp)*Wm-ones(size(Wp))) < 1e-15
        continue;
    end
    [~,s,t]=CheckEnv(Wp,Wm,fp,neww);
    [~,su]=CheckEnv(Wpu,Wmu,fp,neww);
    tmin(i)=t(find(s>su,1,'last'));
    %generate trial
        qp=qu;
        qm=qu;
        qp(len/2)=pert(i);
        qm(len/2)=pert(i);
        [Wp,Wm,neww]=MakeMultistate(qp,qm);
        [Wpu,Wmu]=MakeMultistate(qu,qu);
    %test trial
    if rcond(fp*Wp+(1-fp)*Wm-ones(size(Wp))) < 1e-15
        continue;
    end
    [~,s,t]=CheckEnv(Wp,Wm,fp,neww);
    [~,su]=CheckEnv(Wpu,Wmu,fp,neww);
    tmax(i)=t(find(s>su,1,'first'));
    if mod(i,10)==0 %|| numulti
        disp(int2str(i));
%         disp(num2str(qp(len/2)));
    end%if
end%for i




end

