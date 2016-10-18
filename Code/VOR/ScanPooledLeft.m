function [ comps ] = ScanPooledLeft( ranges, n )
%comps=SCANPOOLEDLEFT(ranges,n) parameter scan for pooled resource model
%   comps  = learning rate differences: WT_nopre - WT_pre
%   ranges = range for each parameter scan
%   n      = number of states
%   parametrs: Pot, Dep_min, 

builder_h = @PooledBuilder;

vexpt=VORbuilder(builder_h, n, ranges(1)*[1 1], ranges(1:2), ranges(1:2), ranges(2), ranges(1), ranges(3), 1,1, true);
range_ctr = 1:length(ranges);

comps = -Inf(length(ranges) * ones(1,6));

for i1 = range_ctr
    Wp = builder_h(n, ranges(i1));
    vexpt.WT = vexpt.WT.setWp(Wp);
    
    for i2 = range_ctr
        for i3 = i2+1:range_ctr(end)
            [~,Wm] = builder_h(n, ranges([i2 i3]));
            vexpt.WT = vexpt.WT.setWm(Wm);
            
            for i4 = range_ctr
                vexpt.nopre = vexpt.nopre.setFp(ranges(i4),1);
                
                for i5 = 1:i4-1
                    vexpt.nopre = vexpt.nopre.setFp(ranges(i5),2);
                    vexpt.withpre = vexpt.withpre.setFp(ranges(i5),3);
                    
                    for i6 = i4+1:range_ctr(end)
                        vexpt.withpre = vexpt.withpre.setFp(ranges(i6),2);
                        
                        comps(i1,i2,i3,i4,i5,i6) = vexpt.InitRComp_left();
                    end%for i6
                end%for i6
            end%for i6
        end%for i6
    end%for i6
end%for i6


end

