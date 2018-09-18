function [ dwelltimes ] = DwellTimes( obj,varargin )
%dwelltimes=obj.DWELLTIMES times dwelt in set of states with same synaptic
%weight
%   dwelltimes = cell, one for each value of synaptic weight, with row
%   vector of dwell times.
%   obj = SynapsePlastSeq, if nonscalar we combine across them


persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='SynapsePlastSeq.DwellTimes';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addOptional('MinNumWvals',0,@(x)validateattributes(x,{'numeric'},{'scalar','integer'},'SynapsePlastSeq.DwellTimes','MinNumWvals',2));
end
p.parse(varargin{:});


dwelltimes=cell(1,max(obj.NumWvals,p.Results.MinNumWvals));

if isscalar(obj)
    wchange=find(diff(obj.readouts)~=0);
    wbefore=obj.readouts(wchange);
    dwell=diff([0 wchange]);
    for i=1:length(dwelltimes)
        dwelltimes{i}=dwell(wbefore==i);
    end%for i
else
    for i=1:numel(obj)
        tempdwell=obj(i).DwellTimes(length(dwelltimes));
        dwelltimes=cellfun(@(x,y) [x y],dwelltimes,tempdwell,'UniformOutput',false);
    end%for i
end%if scalar


end

