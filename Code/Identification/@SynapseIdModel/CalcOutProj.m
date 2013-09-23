function [ obj ] = CalcOutProj( obj )
%obj=CALCOUTPROJ(obj) calculate obj.outProj from obj.w
%   w       = synaptic weights.
%   outProj = cell of diagonal matrices for each possible value of
%             output(low to high), with elements equal to prob of output

wvals=obj.GetWVals;

for i=1:length(wvals)
    obj.outProj{i}=diag(obj.w==wvals(i));
end

obj.outProj(i+1:end)=[];

end

