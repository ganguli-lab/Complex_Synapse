classdef SynapsePlastSeq
    %SYNAPSEPLASTSEQ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess=protected)%data
        %seq of plast types (indices of SynapseIdModel.M). default: []
        potdep=[],
        %sequence of synaptic states. default: [1]
        stateseq=1;
        %sequence synaptic weights (indices of SynapseIdModel.outProj). default: [1]
        readouts=1;
    end
    
    methods (Access=?SynapseIdModel) %setting data
        %
        function newobj=setPotDep(obj,newPotdep)
            newobj=obj;
            newobj.potdep=newPotdep;
        end
        %
        function newobj=setStateSeq(obj,newStateSeq)
            newobj=obj;
            newobj.stateseq=newStateSeq;
        end            
        %
        function newobj=setReadouts(obj,newReadouts)
            newobj=obj;
            newobj.readouts=newReadouts;
        end            
    end
    
    methods%validity etc.
        tf=isvalid(obj)
        tf=iscompatible(obj,modelobj)
    end
    
    methods
        newobj=GetRange(obj,range)
    end
    
    methods (Access=private)%for constructiuon
        %called by constructor
        copy=CopyProps(original,copy)
        [s,x] = assignToObject(s, x)
    end%methods
    
    methods%constructor
        function obj=SynapsePlastSeq(varargin)
            if nargin ~=0%false -> default constructor does nothing
                if nargin==2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
                    %true -> preallocate with default constructor doing nothing
                    obj(max(varargin{1},1),max(varargin{2},1))=SynapsePlastSeq;
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
                    [tempobj,varargin]=extractArgOfType(varargin,'SynapsePlastSeq');
                    %otherwise
                    if isempty(tempobj)
                        tempobj=obj;
                    end
                    %
                    %Set size of object:
                    %
                    if nargin>=2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
                        obj(varargin{1},varargin{2})=SynapsePlastSeq;
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
        end%function SynapsePlastSeq
    end%methods constructor
    
end

