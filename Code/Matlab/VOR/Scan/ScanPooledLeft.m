function [ comps ] = ScanPooledLeft( ranges, n )
%comps=SCANPOOLEDLEFT(ranges,n) parameter scan for pooled resource model
%   comps  = learning rate differences: WT_nopre - WT_pre
%   ranges = range for each parameter scan
%   n      = number of states
%   parametrs: pot, dep_min, dep_max, fp_norm, fp_inc, fp_dec

builder_h = @PooledBuilder;

vexpt=VORbuilderKO(builder_h, n, ranges(1)*[1 1], ranges(1:2), ranges(1:2), ranges(2), ranges(1), ranges(3), 1,1, true);
m = length(ranges);

comps = -Inf(length(ranges) * ones(1,6));

for i1 = 1:m
    DispCounter(i1,m,'i1:');
    Wp = builder_h(n, ranges(i1));
    vexpt.WT = vexpt.WT.setWp(Wp);
    
    DispCounter(1,m,'i2:');
    for i2 = 2:m
        DispCounter(i2,m,'i2:');
        for i3 = 1:i2-1
%             DispCounter(i3,i2-1,'i3:');
            [~,Wm] = builder_h(n, ranges([i3 i2]));
            vexpt.WT = vexpt.WT.setWm(Wm);
            
%             DispCounter(1,m-1,'i4:');
            for i4 = 2:m-1
%                 DispCounter(i4,m-1,'i4:');
                vexpt.nopre = vexpt.nopre.setFp(ranges(i4),1);
                
                for i5 = 1:i4-1
%                     DispCounter(i1,i4-1,'i5:');
                    vexpt.nopre = vexpt.nopre.setFp(ranges(i5),2);
                    vexpt.withpre = vexpt.withpre.setFp(ranges(i5),3);
                    
%                     DispCounter(1,m,'i6:');
                    for i6 = i4+1:m
%                         DispCounter(i6,m,'i6:');
                        vexpt.withpre = vexpt.withpre.setFp(ranges(i6),2);
                        
                        comps(i1,i2,i3,i4,i5,i6) = vexpt.InitRComp_left();
                    end%for i6
%                     DispCounter(m+1,m,'i6:');
                end%for i5
%                 DispCounter(i4,i4-1,'i5:');
            end%for i4
%             DispCounter(m,m-1,'i4:');
        end%for i3
%         DispCounter(i2,i2-1,'i3:');
    end%for i2
    DispCounter(m+1,m,'i2:');
end%for i2
DispCounter(m+1,m,'i1:');


end

