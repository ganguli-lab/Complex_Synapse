function [ comps, wh ] = ScanCNtop2( ranges, n, useCnotN, varargin  )
%comps=SCANCNTOP(ranges,n,useCnotN) parameter scan for cascade/nonuni
%multistate model
%   comps  = learning rate gradient wrt qdep
%   ranges = range for each parameter scan
%   n      = number of states
%   useCnotN = true for cascade, false for nonuniform multistate model
%   parametrs: pot_wt, dep_wt, fp_norm,

if useCnotN
    builder_h = @CascadeBuilder;
    grad_h = @CascadeGrad;
    fr = 2;
else
    builder_h = @NonuniBuilder;
    grad_h = @NonuniGrad;
    fr = 1;
end

modelobj = SynapseMemoryModel.Build(builder_h, ranges(1), {n, ranges(1)});
range_ctr = 1:length(ranges);
m = range_ctr(end);

comps = NaN(m * ones(1,3));
wh=[];

for i1 = range_ctr
%     DispCounter(i1,m,'i1:');
    
    Wp = builder_h(n, ranges(i1));
    dWp = grad_h(n, ranges(i1));
    modelobj = modelobj.setWp(Wp);
    
    for i2 = range_ctr
%         DispCounter(i2,m-1,'i2:');
        
        [~,Wm] = builder_h(n, ranges(i2));
        [~,dWm] = grad_h(n, ranges(i2));
        modelobj = modelobj.setWm(Wm);
            
        for i3 = range_ctr
%             DispCounter(i3,m,'i3:');

            modelobj = modelobj.setFp(ranges(i3)*fr);

            comps(i1,i2,i3) = HorizDeriv(modelobj, dWp, dWm, true);
            if comps(i1,i2,i3) < 0
                wh = [wh; i1, i2, i3, comps(i1,i2,i3)];
            end
 
        end%for i3
%         DispCounter(m+1,m,'i3:');
    end%for i2
%     DispCounter(m+1,m,'i2:');
end%for i1
% DispCounter(m+1,m,'i1:');

end

