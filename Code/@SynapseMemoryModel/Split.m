function [ objs,fracs ] = Split( obj,ix )
%[objs,fracs]=SPLIT(obj,ix) Split SynapseMemoryModel obj into objs
%   ix = logical indexes for objs(1)
%   objs = vector of SynapseMemoryModels
%   fracs = probability of being in each model


objs(1)=SynapseMemoryModel('Wp',obj.Wp(ix,ix),'Wm',obj.Wm(ix,ix),'fp',obj.fp,'w',obj.w(ix));
objs(2)=SynapseMemoryModel('Wp',obj.Wp(~ix,~ix),'Wm',obj.Wm(~ix,~ix),'fp',obj.fp,'w',obj.w(~ix));

p1=objs(1).EqProb;
p2=objs(3).EqProb;

f12=sum(p1*W(ix,~ix));
f21=sum(p2*W(~ix,ix));

fracs(1)=f21/(f12+f21);
fracs(2)=f12/(f12+f21);

end

