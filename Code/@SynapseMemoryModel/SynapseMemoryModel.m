classdef SynapseMemoryModel
    %SYNAPSEMEMORYMODEL Model of complex synapse. For computing memory
    %curves, etc.
    %   Detailed explanation goes here
    
    properties (SetAccess=protected)%data
        %cts time Markov matrix for potentiation
        Wp=zeros(2);
        %cts time Markov matrix for depression
        Wm=zeros(2);
        %synaptic weights. default: [-1; 1]
        w=[-1; 1];
        %fraction of candidate plasticity events that are potentiating
        fp=0.5;
    end
    
    properties%thresholds
        %smallest singular value for Split models
        SingValThresh=1e-10;
        %smallest reciprocal condition number for inverting Zinv
        RCondThresh=1e-5;
        %degeneracy threshold, for evals or eta^+
        DegThresh=1e-3;
        %threshold for lumpability test
        LumpThresh=1e-3;
        %threshold for orphaned states
        OrphanThresh=1e-3;
    end
    
    methods %setting data
        %
        function newobj=setWp(obj,newWp)
            newobj=obj;
            newobj.Wp=newWp;
        end            
        %
        function newobj=setWm(obj,newWm)
            newobj=obj;
            newobj.Wm=newWm;
        end
        %
        function newobj=setW(obj,newW)
            newobj=obj;
            newobj.w=newW;
        end            
        %
        function newobj=setFp(obj,newFp)
            newobj=obj;
            newobj.fp=newFp;
        end
    end
    
    methods%validity etc.
        newobj=Normalise(obj)
        newobj=Zero(obj)
        [newobj,ix]=Sort(obj)
        newobj=Reorder(obj,ix)
        tf=isvalid(obj)
    end
    
    methods %size info
        %
        function val=NumStates(obj)
            val=length(obj.w);
        end
        %
        tf=SameSizes(obj,otherobj)
    end
    
    methods%calculations etc.
        obj3=plus(obj1,obj2)
        obj3=minus(obj1,obj2)
        obj3=mtimes(obj1,obj2)
        obj3=mrdivide(obj1,obj2)
    end
    
    methods%associated matrices
        Wf=GetWf(obj)
        q=GetEnc(obj)
        [Zinv,piv]=GetZinv(obj,varargin)
        p=EqProb(obj,varargin)
        mw=MeanWeight(obj,varargin)
    end
    
    methods%for memory curves etc
        [qa,ca]=Spectrum(obj,varargin)
        S=SNRcurve(obj,t,varargin)
        snr=SNRlearn(obj,t,rlearn,tlearn)
        area=SNRarea(obj)
        A=SNRlaplace(obj,s)
        A=SNRrunAve(obj,t)
        T=FPT(obj)
        deta=DeltaEta(obj)
        tau=MixTime(obj)
        [As,dAp,dAm]=SNRlaplaceGrad(obj,s)
    end
    
    methods%simplifying
        [tf,wtest,Mtest]=TestLump(obj,partitions,varargin)
        partitions=FindLumps(obj,varargin)
        newobj=Lumpify(obj,partitions)
    end
    
    methods (Static=true) %for construction
        newobj=Build(func,fp,varargin)
        newobj=Rand(w,fp,varargin)
        memmodelobj=FromIdModel(idmodelobj)
    end
    
    methods (Access=private)%for constructiuon
        %called by constructor
        copy=CopyProps(original,copy)
        [s,x] = assignToObject(s, x)
    end%methods
    
    methods%constructor
        function obj=SynapseMemoryModel(varargin)
            if nargin ~=0%false -> default constructor does nothing
                if nargin==2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
                    %true -> preallocate with default constructor doing nothing
                    obj(max(varargin{1},1),max(varargin{2},1))=SynapseMemoryModel;
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
                    [tempobj,varargin]=extractArgOfType(varargin,'SynapseMemoryModel');
                    %otherwise
                    if isempty(tempobj)
                        tempobj=obj;
                    end
                    %
                    %Set size of object:
                    %
                    if nargin>=2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
                        obj(varargin{1},varargin{2})=SynapseMemoryModel;
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


