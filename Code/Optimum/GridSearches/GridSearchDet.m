function [ results_struct ] = GridSearchDet( nrange, trange,  varargin )
%GRIDSEARCH Summary of this function goes here
%   Detailed explanation goes here



results_struct(length(nrange),length(trange))=struct('Wp',[],'Wm',[],'SNR',[]);


for a=1:length(nrange)
    n=nrange(a);
    w=BinaryWeights(n);
%     targets=zeros(numtrials,n-1);
%     for i=1:n-1
%             targets(:,i)=i+randi(n-i,1,numtrials);
%     end%for i
    for b=1:length(trange)
        smax=0;
        scurve=zeros(1,length(trange));
        for i1=1:n-1
            for j1=i1+1:n
        for i2=i1+1:n-1
            for j2=i2+1:n
        for i3=i2+1:n-1
            for j3=i3+1:n
%         for i4=i3+1:n-1
%             for j4=i4+1:n
                targets=[2:n n];
                targets(i1)=j1;
                targets(i2)=j2;
                targets(i3)=j3;
%                 targets(i4)=j4;
                [nWp,nWm]=Deterministic(targets);
                snew=SNRcurve(trange,nWp,nWm,0.5,w,'UseExpM',true);
                if any(snew>1)
                    continue;
                end
                if snew(b)>smax
                    Wp=nWp;
                    Wm=nWm;
                    smax=snew(b);
                    scurve=snew;
                end%if
%             end%for j4
%         end%for i4
            end%for j3
        end%for i3
            end%for j2
        end%for i2
            end%for ji
        end%for i1
        if n==2
            [Wp,Wm]=Deterministic([2 2]);
            scurve=SNRcurve(trange,Wp,Wm,0.5,w,'UseExpM',true);
        end
        results_struct(a,b)=struct('Wp',Wp,'Wm',Wm,'SNR',scurve);
    end%for b
    disp(n);
end%for a



end

