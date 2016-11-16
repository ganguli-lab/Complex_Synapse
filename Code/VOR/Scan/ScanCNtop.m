function [ comps ] = ScanCNtop( ranges, n, reps, useCnotN  )
%comps=SCANCNTOP(ranges,n,reps,useCnotN) parameter scan for cascade/nonuni
%multistate model
%   comps  = learning rate differences: WT_nopre - KO_nopre
%   ranges = range for each parameter scan
%   n      = number of states
%   useCnotN = true for cascade, false for nonuniform multistate model
%   parametrs: pot_wt, dep_wt, dep_ko, fp_norm, fp_inc

if useCnotN
    builder_h = @CascadeBuilder;
    maxv = 0.5;
else
    builder_h = @NonuniBuilder;
    maxv = 1;
end

vexpt=VORbuilderKO(builder_h, n, ranges(1), ranges(1), ranges(2), ranges(2), ranges(1), 0.5, 1,1, false);
range_ctr = 1:length(ranges);
m = range_ctr(end);

comps = -Inf(length(ranges) * ones(1,5));

for i1 = range_ctr
    DispCounter(i1,m,'i1:');
    
    Wp = builder_h(n, ranges(i1));
    vexpt.WT = vexpt.WT.setWp(Wp);
    
    for i2 = range_ctr(1:end-1)
        DispCounter(i2,m-1,'i2:');
        
        [~,Wm] = builder_h(n, ranges(i2));
        vexpt.WT = vexpt.WT.setWm(Wm);
        
        DispCounter(1,m,'i3:');
        for i3 = range_ctr(range_ctr > i2)
%             pause;
            DispCounter(i3,m,'i3:');
%             pause;
            
            [~,Wm] = builder_h(n, ranges(i3));
            vexpt.KO = vexpt.KO.setWm(Wm);
            
            DispCounter(1,m,'i4:');
            for i4 = range_ctr(2:end)
                DispCounter(i4,m,'i4:');
                
                vexpt.nopre = vexpt.nopre.setFp(ranges(i4)/maxv,1);
                
                x = Find_pot_KO(builder_h,n,ranges(i1),ranges(i2),ranges(i3),ranges(i4)/maxv,reps,0,maxv);
                Wp = builder_h(n, x);
                vexpt.KO = vexpt.KO.setWp(Wp);
                
                for i5 = range_ctr(range_ctr < i4)
%                     DispCounter(i5,i4-1,'i5:');
                    
                    vexpt.nopre = vexpt.nopre.setFp(ranges(i5)/maxv,2);

                    comps(i1,i2,i3,i4,i5) = vexpt.InitRComp_top();
%                       pause;
                end%for i5
%                 DispCounter(i4,i4-1,'i5:');
            end%for i4
%             pause;
            DispCounter(m+1,m,'i4:');
%             pause;
        end%for i3
%          pause;
        DispCounter(m+1,m,'i3:');
%          pause;
    end%for i2
    DispCounter(m,m-1,'i2:');
end%for i1
DispCounter(m+1,m,'i1:');

end

