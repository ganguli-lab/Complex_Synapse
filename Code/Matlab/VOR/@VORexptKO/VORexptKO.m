classdef VORexptKO < VORexperiment
    %VOREXPTKO class for plotting VOR comparisons between WT/KO
    %with/without pretraining
    %   Detailed explanation goes here
    
    properties
        %SynapseMemoryModel for MHC-I DbKb knockout
        KO=SynapseMemoryModel;
    end
    
    properties %labels
        KOlabel='DKO';
        KOcolor=[192 0 0]/255;
        KOstyle='-';
        %
    end
    
    methods
        ProbEvols( obj,fh,varargin )
        EqProbPlots( obj,fh,varargin )
        PlotLearn( obj,varargin )
        PlotLearnS( obj,varargin )
        PrintFigs( obj,prefix )
        St=LearnSdata(obj,varargin);
        [P_WT_nopre,P_KO_nopre,P_WT_pre,P_KO_pre,t]=ProbEvolsData(obj)
        comps=InitialRateComps(obj)
    end
    
    methods (Access=private)%for constructiuon
        %called by constructor
        copy=CopyProps(original,copy)
        [s,x] = assignToObject(s, x)
    end%methods
    
    methods%constructor
        function obj=VORexptKO(varargin)
            superargs=varargin;
                if nargin>=2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
                    superargs(1:2) = [];
                end%if nargin>=2
            obj@VORexperiment(superargs{:});
            if nargin ~=0%false -> default constructor does nothing
                if nargin==2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
                    %true -> preallocate with default constructor doing nothing
                    obj(max(varargin{1},1),max(varargin{2},1))=VORexptKO;
                    if varargin{1}<1
                        obj(1,:)=[];
                    end
                    if varargin{2}<1
                        obj(:,1)=[];
                    end
                else
                    %
                    %default parameters:
                    %if we're copying another obj
                    [tempobj,varargin]=extractArgOfType(varargin,'VORexptKO');
                    %otherwise
                    if isempty(tempobj)
                        tempobj=obj;
                    end
                    %
                    %Set size of object:
                    %
                    if nargin>=2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
                        obj(varargin{1},varargin{2})=VORexptKO;
                    end%if nargin>=2
                    %
                    %set parameter values:
                    [tempobj,varargin]=assignToObject(tempobj,varargin);
                    obj=CopyProps(tempobj,obj);
                    %
                    %Extract data:
                    %
                end% if nargin=2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
            end%if nargin ~=0
        end%function SynapseIdModel
    end%methods constructor
    
end

