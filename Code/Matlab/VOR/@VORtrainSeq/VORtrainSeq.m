classdef VORtrainSeq
    %VORTRAINSEQ class for plotting VOR learning curves
    %   Detailed explanation goes here
    
    properties (SetAccess=protected)%data
        %times at which training regimes end
        tTrain=[];
        %fraction of potentiating events in each training epoch, incl before 
        fps=0.5;
        %rate of candidate plasticity events in each epoch, incl before
        rs=1;
    end
    
    methods %setting data
        %
        function newobj=setT(obj,newT)
            newobj=obj;
            newobj.tTrain=newT;
        end            
        %
        function newobj=setFp(obj,newFp,varargin)
            newobj=obj;
            if isempty(varargin)
                newobj.fps=newFp;
            else
                newobj.fps(varargin{1})=newFp;
            end
        end
    end
    
    methods
        function val=numTrain(obj)
        %number of training epochs
            val=length(obj.tTrain);
        end
        [tf]=isvalid(obj)
        [S,Pt,t]=LearningCurve(obj,modelobj,dt)
        [S,Pt,t]=LearningCurveEnd(obj,modelobj,dt)
        rate=InitialLearnRate(obj,modelobj)
    end
    
    methods (Access=private)%for constructiuon
        %called by constructor
        copy=CopyProps(original,copy)
        [s,x] = assignToObject(s, x)
    end%methods
    
    methods%constructor
        function obj=VORtrainSeq(varargin)
            if nargin ~=0%false -> default constructor does nothing
                if nargin==2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
                    %true -> preallocate with default constructor doing nothing
                    obj(max(varargin{1},1),max(varargin{2},1))=VORtrainSeq;
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
                    [tempobj,varargin]=extractArgOfType(varargin,'VORtrainSeq');
                    %otherwise
                    if isempty(tempobj)
                        tempobj=obj;
                    end
                    %
                    %Set size of object:
                    %
                    if nargin>=2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
                        obj(varargin{1},varargin{2})=VORtrainSeq;
                    end%if nargin>=2
                    %
                    %set parameter values:
                    [tempobj,varargin]=assignToObject(tempobj,varargin);
                    obj=CopyProps(tempobj,obj);
                    %
                    %Extract data:
                    %
                end% if nargin=2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
                if length(obj.rs) <= length(obj.fps)
                    obj.rs = [obj.rs ones(1, length(obj.fps) - length(obj.rs))];
                end%if length(rs)
            end%if nargin ~=0
        end%function SynapseIdModel
    end%methods constructor
    
end

