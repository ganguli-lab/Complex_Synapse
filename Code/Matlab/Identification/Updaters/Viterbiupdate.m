function [ newmodelobj,loglike,pstate ] = Viterbiupdate( modelobj,simobj,varargin )
%[newmodelobj,loglike,pstate]=VITERBIUPDATE(modelobj,simobj) Viterbi update of estimated HMM
%   newmodelobj = updated SynapseIdModel
%   loglike     = log likelihood of readout given current model for most likely path
%   pstate      = posterior prob of HMM being in each state at each time, delta fn at most likely path
%   modelobj = SynapseIdModel
%   simobj   = SynapsePlastSeq



persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='Viterbiupdate';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addParameter('Normalise',true);
%     p.addParameter('Normalise',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
end
p.parse(varargin{:});


M_new={zeros(length(modelobj.M{1}))};
for i=2:numel(modelobj.M)
    M_new{i}=M_new{1};
end

loglike=log(modelobj.Initial*modelobj.outProj{simobj.readouts(1)});

LikelyPath=zeros(modelobj.NumStates,simobj.NumT);
ev=ones(1,modelobj.NumStates);
n=(1:modelobj.NumStates)';
pstate=LikelyPath;

M=cat(3,modelobj.M{simobj.potdep(1:end-1)});
outProj=cat(3,modelobj.outProj{simobj.readouts(2:end)});

updater=log(mmx('mult',M,outProj));


for i=1:simobj.NumT-1
    ExpandPath(i);
end

LikelyPath(:,end)=n;
[loglike,ix]=max(loglike);
LikelyPath=LikelyPath(ix,:);

stepstart=LikelyPath(1:end-1);
stepend=LikelyPath(2:end);
steps=sub2ind(modelobj.NumStates*[1 1],stepstart,stepend);

for i=1:modelobj.NumPlast
    M_new{i}(1:modelobj.NumStates^2)=hist(steps(simobj.potdep(1:end-1)==i),1:modelobj.NumStates^2);
end
    
% for t=1:simobj.NumT-1
%     M_new{simobj.potdep(t)}(LikelyPath(t),LikelyPath(t+1))=M_new{simobj.potdep(t)}(LikelyPath(t),LikelyPath(t+1)) + 1;
% end

for i=1:length(M_new)
    M_new{i}=M_new{i}+diag(sum(M_new{i},2)==0);
end

pstate(sub2ind(size(pstate),LikelyPath,1:simobj.NumT))=1;

newmodelobj=modelobj.setM(M_new);
% newmodelobj=newmodelobj.setInitial(pstate(:,1)');
newmodelobj=newmodelobj.setInitial(sum(pstate,2)');

if p.Results.Normalise
    newmodelobj=newmodelobj.Normalise;
    assert(newmodelobj.isvalid,'newmodelobj is invalid');
end




    function ExpandPath(step)
        LikelyPath(:,step)=n;
        StepLike = loglike'*ev + updater(:,:,step);
        [loglike,ix]=max(StepLike,[],1);
        LikelyPath=LikelyPath(ix,:);
    end

        

end

