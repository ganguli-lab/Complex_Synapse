function [ results_struct ] = GridSearchShortcut( nrange, trange, qsrange, varargin )
%GRIDSEARCH Summary of this function goes here
%   Detailed explanation goes here



results_struct(length(nrange),length(trange))=struct('Wp',[],'Wm',[],'SNR',[]);


for a=1:length(nrange)
    n=nrange(a);
    w=BinaryWeights(n);
    for b=1:length(trange)
        smax=0;
        scurve=zeros(1,length(trange));
        for qs=qsrange
             for i=1:(n-1)
                for j=(i+1):n
                    for k=1:length(qsrange)
                        qd=qsrange(k)*(1-qs);
                        [nWp,nWm]=OneShortcut(n,i,j,qs,qd);
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
                    end%for k
                end%for j
            end%for i
        end%for qs
        results_struct(a,b)=struct('Wp',Wp,'Wm',Wm,'SNR',scurve);
    end%for b
    disp(n);
end%for a



end

