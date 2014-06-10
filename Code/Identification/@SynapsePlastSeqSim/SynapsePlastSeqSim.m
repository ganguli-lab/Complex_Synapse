classdef SynapsePlastSeqSim < SynapsePlastSeq
    %SYNAPSEPLASTSEQSIM Simulated sequence of plasticity events.
    %   This class stores the results of a simulated sequence of plasticity
    %   events: the type of each event (SYNAPSEPLASTSEQ.potdep), the
    %   synapticstate before/after each event (SYNAPSEPLASTSEQ.stateseq)
    %   and the corresponding synaptic weight index (SYNAPSEPLASTSEQ.readouts)
    
    properties (SetAccess=protected)%data
        %sequence of synaptic states. default: [1]
        stateseq=[];
    end
    
    methods (Access={?SynapseIdModel}) %setting data
        function newobj=setStateSeq(obj,newStateSeq)
            newobj=obj;
            newobj.stateseq=newStateSeq;
        end            
    end
    
    methods%validity etc.
        tf=isvalid(obj)
        tf=iscompatible(obj,modelobj)
    end
    
    methods %size info
        %
        function val=NumStates(obj)
            val=max([obj.stateseq]);
        end
        %
        tf=SameSizes(obj,otherobj)
    end
    
    methods
        newobj=GetRange(obj,range)
        obj3=plus(obj1,obj2)
    end
    
    methods (Access=private)%for constructiuon
        %called by constructor
        copy=CopyProps(original,copy)
        [s,x] = assignToObject(s, x)
    end%methods
    
    methods%constructor
        function obj=SynapsePlastSeqSim(varargin)
            obj@SynapsePlastSeq(varargin{:});
            if nargin ~=0%false -> default constructor does nothing
                if nargin==2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
                    %true -> preallocate with default constructor doing nothing
                    obj(max(varargin{1},1),max(varargin{2},1))=SynapsePlastSeqSim;
                    if varargin{1}<1
                        obj(1,:)=[];
                    end
                    if varargin{2}<1
                        obj(:,1)=[];
                    end
                elseif nargin==1 && isnumeric(varargin{1}) && isrow(varargin{1})
                    siz=num2cell(varargin{1});
                    obj(siz{:})=SynapsePlastSeqSim;
                else
                    %
                    %default parameters:
                    %if we're copying another obj
                    [tempobj,varargin]=extractArgOfType(varargin,'SynapsePlastSeqSim');
                    %otherwise
                    if isempty(tempobj)
                        tempobj=obj;
                    end
                    %
                    %Set size of object:
                    %
                    if nargin>=2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
                        obj(varargin{1},varargin{2})=SynapsePlastSeqSim;
                    elseif nargin>=1 && isnumeric(varargin{1}) && isrow(varargin{1})
                        siz=num2cell(varargin{1});
                        obj(siz{:})=SynapsePlastSeqSim;
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
        end%function SynapsePlastSeqSim
    end%methods constructor
    
end

