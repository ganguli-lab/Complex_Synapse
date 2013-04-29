function [ results_struct ] = GridSearch( nrange, trange, num_trials, varargin )
%GRIDSEARCH Summary of this function goes here
%   Detailed explanation goes here



results_struct(length(nrange),length(trange))=struct('Wp',[],'Wm',[],'SNR',[]);


for i=1:length(nrange)
    n=nrange(i);
    w=BinaryWeights(n);
    for j=1:length(trange)
        t=trange(j);
        s=0;
        for k=1:num_trials
            [ nWp,nWm ] = FindOpt(t,n,varargin{:});
%             [u,d]=eig(0.5*(nWp+nWm));
%             if rcond(u)<1e-5
%                 continue;
%             end
            snew=SNRcurve(trange,nWp,nWm,0.5,w,'UseExpM',true);
            if any(snew>1)
                continue;
            end
            if snew(j)>s
                Wp=nWp;
                Wm=nWm;
                s=snew(j);
            end                
        end
        for k=1:num_trials
            [ nWp,nWm ] = FindOpt(t,n,'Triangular',true,varargin{:});
%             [u,d]=eig(0.5*(nWp+nWm));
%             if rcond(u)<1e-5
%                 continue;
%             end
            snew=SNRcurve(trange,nWp,nWm,0.5,w,'UseExpM',true);
            if any(snew>1)
                continue;
            end
            if snew(j)>s
                Wp=nWp;
                Wm=nWm;
                s=snew(j);
            end                
        end
        results_struct(i,j)=struct('Wp',Wp,'Wm',Wm,'SNR',snew);
    end
    disp(n);
end



end

