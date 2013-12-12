classdef VORexperiment
    %VOREXPERIMENT class for plotting VOR comparisons between WT/KO
    %with/without pretraining
    %   Detailed explanation goes here
    
    properties
        %SynapseMemoryModel for Wild-Type
        WT=SynapseMemoryModel;
        %SynapseMemoryModel for MHC-I DbKb knockout
        KO=SynapseMemoryModel;
        %VORtrainSeq without pre-training
        nopre=VORtrainSeq;
        %VORtrainSeq with pre-training
        withpre=VORtrainSeq;
    end
    
    properties %labels
        WTlabel='WT';
        KOLabel='D^bK^{b-/-}';
        noprelabel='No Pre-training';
        withprelabel='w/ Pre-training';
        WTcolor='k';
        KOcolor=[192 0 0]/255;
        noprestyle='-';
        withprestyle='--';
        %
        LabFontSize=20;
        FontSize=20;
        txFontSize=10;
        EqFontSize=20;
        ProbFontSize=10;
    end
    
    methods
        ProbEvols( obj,fh,varargin )
        EqProbPlots( obj,fh,varargin )
        PlotLearn( obj,varargin )
        PlotLearnS( obj,varargin )
        ViewFigs( obj )
        PrintFigs( obj,prefix )
    end
    
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

