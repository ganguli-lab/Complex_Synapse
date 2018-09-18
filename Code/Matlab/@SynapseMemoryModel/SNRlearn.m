function [ snr ] = SNRlearn( obj, t, rlearn, tlearn )
%snr=SynapseMemoryModel.SNRLEARN(t,rlearn,tlearn) Learning and subsequent
%forgetting curve
%   snr    = SNR(t)
%   t      = vector of time values
%   rlearn = rate of learning events / rate of forgetting events
%   tlearn = time of last learning event

snr=zeros(size(t));

pr = obj.EqProb;
W = obj.GetWf;
Id=eye(size(W));
q = obj.fp*(expm((W + rlearn * obj.Wp) * tlearn) - Id) - (1-obj.fp)*(expm((W + rlearn * obj.Wm) * tlearn) - Id);


ilearn = find(t <= tlearn);
iforget = find(t > tlearn);

for i=ilearn
    snr(i) = pr * (obj.fp*(expm((W + rlearn * obj.Wp) * t(i)) - Id) - (1-obj.fp) * (expm((W + rlearn * obj.Wm) * t(i)) - Id)) * obj.w;
end
for i=iforget
    snr(i) = pr * q * expm(W*(t(i)-tlearn)) * obj.w;
end



end

