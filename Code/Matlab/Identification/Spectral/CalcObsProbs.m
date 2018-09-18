function [ P1,P2,P3 ] = CalcObsProbs( simobj )
%CALCOBSPROBS(simobj) Summary of this function goes here
%   Detailed explanation goes here


P1=hist([simobj.readouts], 1:simobj.NumWvals);
P1=P1/sum(P1);

emptymat=zeros(simobj.NumWvals);
emptyarr=zeros(simobj.NumWvals*[1 1 1]);

P2=cell(1,simobj.NumPlast);
P3=cell(simobj.NumPlast,simobj.NumPlast);

for a=1:length(P2)
    P2{a}=emptymat;
    for b=1:length(P3)
        P3{a,b}=emptyarr;
    end
end

steps=[];
steptypes=[];

dsteps=[];
dsteptypes=[];

for i=1:numel(simobj)

    stepstart=simobj(i).readouts(1:end-1);
    stepend=simobj(i).readouts(2:end);
    steps=[steps sub2ind(size(emptymat),stepstart,stepend)];
    steptypes=[steptypes simobj(i).potdep(1:end-1)];

    stepstart=stepstart(1:end-1);
    stepmid=stepend(1:end-1);
    stepend=stepend(2:end);
    dsteps=[dsteps sub2ind(size(emptyarr),stepstart,stepmid,stepend)];
    dsteptypes=[dsteptypes [simobj(i).potdep(1:end-2); simobj(i).potdep(2:end-1)] ];

end

for a=1:length(P2)
    P2{a}(1:numel(emptymat))=hist(steps(steptypes==a),1:numel(emptymat));
    P2{a}=P2{a}/sum(P2{a}(:));
    for b=1:length(P3)
        P3{a,b}(1:numel(emptyarr))=hist(dsteps(dsteptypes(1,:)==a & dsteptypes(2,:)==b),1:numel(emptyarr));
        P3{a,b}=P3{a,b}/sum(P3{a,b}(:));
    end
end




end

