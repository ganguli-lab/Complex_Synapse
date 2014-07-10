function [ fitmodel ] = FitSynapseSizeDwell( simobjs,varargin )
%[fitmodel,like_n]=FitSynapseSizeDwell(simobjs,synapseOptions,optimOptions)
%Determinig number of states needed for SynapseIdModel (fitmodel) to fit
%SypssePlastSeq (simobjs) using distribution of dwell times in each
%synaptic weight
%   synapseOptions = SynapseOptimset
%   optimOptions  = optimoptions for fmincon
%other arguments passed to SynapseIdModel.Rand

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='FitSynapseSizeDwell';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addOptional('synapseOptions',SynapseOptimset,@(x)validateattributes(x,{'SynapseOptimset'},{},'FitSynapseSizeDwell','synapseOptions',2));
    p.addOptional('optimOptions',optimoptions('fmincon','Algorithm','interior-point'),@(x)validateattributes(x,{'optim.options.Fmincon'},{},'FitSynapseSizeDwell','optimOptions',3));
    p.addOptional('extraArgs',{},@(x)validateattributes(x,{'cell'},{},'FitSynapseSizeDwell','extraArgs',4));
%     p.addParameter('Normalise',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
end
p.parse(varargin{:});
synapseOptions=p.Results.synapseOptions;

% NumPlastTypes=simobjs.NumPlast;
numWvals=simobjs.NumWvals;

dwell=cell(size(simobjs,1),numWvals);
for i=1:size(simobjs,1)
    dwell(i,:)=simobjs(i,:).DwellTimes(numWvals);
end

numStates=zeros(numWvals,1);
for j=1:numWvals
    DispStates(j);
    numStates(j)=FitNumExp( dwell(:,j),synapseOptions,p.Results.optimOptions );
end

numStates=cumsum(numStates)+1;
w=zeros(numStates(end),1);
w(numStates)=1;
w=cumsum(w)+1;
w(end)=[];
fitmodel=SynapseIdModel.Rand(w,p.Results.extraArgs{:});


% like_n=struct('numstates',1:(numWvals*options.MaxStates),'loglike',NaN(1,numWvals*options.MaxStates));
% 
    function DispStates(k)
        if ~strcmpi(synapseOptions.Display,'off')
            disp(['Weight value #' int2str(k)]);
        end
    end


end

