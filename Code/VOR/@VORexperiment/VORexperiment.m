classdef VORexperiment
    %VOREXPERIMENT class for plotting VOR comparisons between WT/KO
    %with/without pretraining
    %   Detailed explanation goes here
    
    properties
        %SynapseMemoryModel for Wild-Type
        WT=SynapseMemoryModel;
        %VORtrainSeq without pre-training
        nopre=VORtrainSeq;
        %VORtrainSeq with pre-training
        withpre=VORtrainSeq;
    end
    
    properties %labels
        WTlabel='WT';
        noprelabel='No Pre-training';
        withprelabel='w/ Pre-training';
        WTcolor='k';
        WTstyle='-';
        noprestyle='-';
        withprestyle='--';
        noprecolor='k';
        withprecolor='c';
        %
        LabFontSize=20;
        FontSize=20;
        LegFontSize=20;
        txFontSize=10;
        EqFontSize=20;
        ProbFontSize=10;
        %
        numpts=100;
        pooled=false;
    end
    
    methods
        ProbEvols( obj,fh,varargin )
        EqProbPlots( obj,fh,varargin )
        PlotLearn( obj,varargin )
        PlotLearnS( obj,varargin )
        ViewFigs( obj )
        PrintFigs( obj,prefix )
        St=LearnSdata(obj,varargin);
        [P_WT_nopre,P_WT_pre,t]=ProbEvolsData(obj)
        comp=InitRComp_left(obj)
    end
    
    methods (Access=private)%for constructiuon
        %called by constructor
        copy=CopyProps(original,copy)
        [s,x] = assignToObject(s, x)
    end%methods
    
    methods%constructor
        function obj=VORexperiment(varargin)
            if nargin ~=0%false -> default constructor does nothing
                if nargin==2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
                    %true -> preallocate with default constructor doing nothing
                    obj(max(varargin{1},1),max(varargin{2},1))=VORexperiment;
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
                    [tempobj,varargin]=extractArgOfType(varargin,'VORexperiment');
                    %otherwise
                    if isempty(tempobj)
                        tempobj=obj;
                    end
                    %
                    %Set size of object:
                    %
                    if nargin>=2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
                        obj(varargin{1},varargin{2})=VORexperiment;
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

