function [ options ] = SynapseOptimset( varargin )
%options=SYNAPSEOPTIMSET(oldoptions,param,val,...) Set optimisation parameters for
%synapse fitting
%   oldoptions can be omitted, defaults will be used instead

options=struct('MaxIter',1000,'TolFun',NaN,'TolX',1e-4,'TolFunChange',1,...
    'Algorithm','BW','Weighter','RJ','ExtraParams',{{}},...
    'Display','off','OutputFcn',[],'PlotFcn',[],...
    'fp',0.5,'GroundTruth',[]);

if isstruct(varargin{1})
    options=UpdateOldOptions(options,varargin{1});
    varargin(1)=[];
end

if nargin>1
    varargin=extractPVpairs(varargin);
    [options,unused]=UpdateOptions(options,varargin{:});
    options.ExtraParams=[options.ExtraParams unused];
end


    function [opt,unused]=UpdateOptions(oldopt,newoptcell)
        opt=oldopt;
        unused=cell(1,length(newoptcell));
        nextun=1;
        for ii=1:2:length(newoptcell)-1
            if isfield(opt,newoptcell{ii})
                opt.(newoptcell{ii})=newoptcell{ii+1};
            else
                unused{nextun}=newoptcell{ii};
                unused{nextun+1}=newoptcell{ii+1};
                nextun=nextun+2;
            end%if isfield
        end%for ii
        unused(nextun:end)=[];
    end

    function [opt,unused]=UpdateOldOptions(oldopt,newopt)
        opt=oldopt;
        pr=fields(newopt);
        unused=cell(1,2*length(pr));
        nextun=1;
        for ii=1:length(pr)
            if ~isempty(newopt.(pr{ii}))
                if isfield(opt,pr{ii})
                    opt.(pr{ii})=newopt.(pr{ii});
                else
                    unused{nextun}=pr{ii};
                    unused{nextun+1}=newopt.(pr{ii});
                    nextun=nextun+2;
                end%if isfield
            end%if ~isempty
        end%for ii
        unused(nextun:end)=[];
    end


end

