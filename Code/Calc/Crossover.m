function [tmin,tmax,pert]=Crossover( varargin )
%CROSSOVER Summary of this function goes here
%   Detailed explanation goes here

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='Crossover';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addParameter('fp',0.5,@(x)validateattributes(x,{'numeric'},{'scalar'},'Crossover','fp'));
    p.addParameter('len',20,@(x)validateattributes(x,{'numeric'},{'scalar','integer'},'Crossover','len'));
    p.addParameter('numtrials',100,@(x)validateattributes(x,{'numeric'},{'scalar','integer'},'Crossover','numtrials'));
    p.addParameter('di',1,@(x)validateattributes(x,{'numeric'},{'scalar'},'Crossover','di'));
    p.addParameter('sh',0.1,@(x)validateattributes(x,{'numeric'},{'scalar'},'Crossover','sh'));
end
p.parse(varargin{:});
r=p.Results;


qu=r.di*ones(1,r.len-1);

pert=zeros(1,r.numtrials);
tmin=pert;
tmax=pert;

for i=1:r.numtrials
    pert(i)=r.sh+(r.di-r.sh)*(r.numtrials-i)/(r.numtrials);
    %generate trial
        qp=qu;
        qm=qu;
        qp(r.len-1)=pert(i);
        qm(1)=pert(i);
        [Wp,Wm,neww]=MakeMultistate(qp,qm);
        [Wpu,Wmu]=MakeMultistate(qu,qu);
    %test trial
    if rcond(r.fp*Wp+(1-r.fp)*Wm-ones(size(Wp))) < 1e-15
        continue;
    end
    [~,s,t]=CheckEnv(Wp,Wm,r.fp,neww);
    [~,su]=CheckEnv(Wpu,Wmu,r.fp,neww);
    tmin(i)=t(find(s>su,1,'last'));
    %generate trial
        qp=qu;
        qm=qu;
        qp(r.len/2)=pert(i);
        qm(r.len/2)=pert(i);
        [Wp,Wm,neww]=MakeMultistate(qp,qm);
        [Wpu,Wmu]=MakeMultistate(qu,qu);
    %test trial
    if rcond(r.fp*Wp+(1-r.fp)*Wm-ones(size(Wp))) < 1e-15
        continue;
    end
    [~,s,t]=CheckEnv(Wp,Wm,r.fp,neww);
    [~,su]=CheckEnv(Wpu,Wmu,r.fp,neww);
    tmax(i)=t(find(s>su,1,'first'));
    if mod(i,10)==0 %|| numulti
        r.disp(int2str(i));
%         r.disp(num2str(qp(r.len/2)));
    end%if
end%for i




end

