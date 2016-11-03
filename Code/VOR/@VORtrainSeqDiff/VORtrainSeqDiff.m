classdef VORtrainSeqDiff
    %VORtrainSeqDiff class for plotting VOR learning curves
    %   Includes compensating synapses
    
    properties (SetAccess=protected)%data
        %Training sequence of relevant synapses for VOR increase
        VORrel = VORtrainSeq;
        %Training sequence of synapses that cancel VOR increase
        VORcomp = VORtrainSeq;
        %fraction of total synapse poopulation in other population of synapses
        frac_other = 0.5;
    end
    
    properties (Dependent=true)
        tTrain;
    end
    
    methods %setting/getting data
        %
        function newobj=setFrac(obj,newF,varargin)
            newobj=obj;
            if isempty(varargin)
                newobj.frac_other=newF;
            else
                newobj.frac_other(varargin{1})=newF;
            end
        end
        %
        function value = get.tTrain(obj)
            value = obj.VORrel.tTrain;
        end
        %
        function obj = set.tTrain(obj, value)
            obj.VORrel = obj.VORrel.setT(value);
            obj.VORcomp = obj.VORcomp.setT(value);
        end
    end
    
    methods
        function val=numTrain(obj)
        %number of training epochs
            val=length(obj.tTrain);
        end
        tf=isvalid(obj)
        [S,Pt,t,Pt_other]=LearningCurve(obj,modelobj,dt)
        [S,Pt,t]=LearningCurveEnd(obj,modelobj,dt)
        rate=InitialLearnRate(obj,modelobj)
    end
    
    methods (Access=private)%for constructiuon
        %called by constructor
        copy=CopyProps(original,copy)
        [s,x] = assignToObject(s, x)
    end%methods
    
    methods%constructor
        function obj=VORtrainSeqDiff(varargin)
%             superargs=varargin;
%                 if nargin>=2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
%                     superargs(1:2) = [];
%                 end%if nargin>=2
%             obj@VORtrainSeq(superargs{:});
            if nargin ~=0%false -> default constructor does nothing
                if nargin==2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
                    %true -> preallocate with default constructor doing nothing
                    obj(max(varargin{1},1),max(varargin{2},1))=VORtrainSeqDiff;
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
                    [tempobj,varargin]=extractArgOfType(varargin,'VORtrainSeqDiff');
                    %otherwise
                    if isempty(tempobj)
                        tempobj=obj;
                    end
                    %
                    %Set size of object:
                    %
                    if nargin>=2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
                        obj(varargin{1},varargin{2})=VORtrainSeqDiff;
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

