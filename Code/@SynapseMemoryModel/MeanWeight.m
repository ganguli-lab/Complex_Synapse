function [ mw ] = MeanWeight( obj,varargin )
%mw=obj.MEANWEIGHT(...) Mean synaptic weight in equilibrium
%   mw  = mean weight
%   obj = SynapseMemoryModel
%   ... passed to obj.EqProb(...)

p = obj.EqProb(varargin{:});
mw = p * obj.w;


end

