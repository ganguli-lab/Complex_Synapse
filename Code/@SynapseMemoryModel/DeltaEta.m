function [ deta ] = DeltaEta( obj )
%deta=SynapseMemoryModel.DELTAETA eta^+_i - eta^+_n
%   W = transition rates
%   w = Weights of states (+/-1)


deta = -obj.GetZinv \ obj.w;

deta=(deta-deta(end))/2;


end

