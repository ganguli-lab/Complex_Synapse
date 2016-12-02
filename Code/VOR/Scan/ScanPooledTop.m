function [ comps,wh ] = ScanPooledTop( prm, n, reps, varargin  )
%comps=SCANCNTOP(ranges,n,reps) parameter scan for pooled
%resource model (dep only)
%   comps  = learning rate differences: WT_nopre - KO_nopre
%   ranges = range for each parameter scan
%   n      = number of states
%   reps   = number of attempts
%   parametrs: pot_wt, dep_wt, dep_ko, fp_norm, fp_inc

builder_h = @PooledBuilder;
grad_h = @PooledGrad;
minv = 1e-4;
maxv = 1 - minv;

vexpt=VORbuilderKO(builder_h, n, prm(1), prm(1), prm(2), prm(2), prm(1), 0.5, 1,1, false);
m = 1:length(prm);
M = m(end);
wh=[];

comps = -Inf(length(prm) * ones(1,7));

for i1 = m
    DispCounter(i1,M,'i1:');
    
    Wp = builder_h(n, prm(i1));
    vexpt.WT = vexpt.WT.setWp(Wp);
    
    for i2 = m(1:end-2)
        DispCounter(i2,M-2,'i2:');
        DispCounter(1,M-1,'i2b:');
        for i2b = m(m > i2 & m < M)
            DispCounter(i2b,M-1,'i2b:');

            [~,Wm] = builder_h(n, [prm(i2b) prm(i2)]);
            vexpt.WT = vexpt.WT.setWm(Wm);

            DispCounter(1,M-1,'i3:');
            for i3 = m(m > i2 & m < M)
                DispCounter(i3,M-1,'i3:');
                DispCounter(1,M,'i3b:');
                for i3b = m(m > i2b & m > i3)
                    DispCounter(i3b,M,'i3b:');

                    [~,Wm] = builder_h(n, [prm(i3b) prm(i3)]);
                    vexpt.KO = vexpt.KO.setWm(Wm);

                    DispCounter(1,M,'i4:');
                    for i4 = m(2:end)
                        DispCounter(i4,M,'i4:');

                        vexpt.nopre = vexpt.nopre.setFp(prm(i4),1);

                        x = Find_pot_KO(builder_h,n,prm(i1),[prm(i2b) prm(i2)],[prm(i3b) prm(i3)],prm(i4),reps,minv,maxv,'ObjGrad',grad_h,varargin{:});
                        if isnan(x)
                            continue;
                        end
                        Wp = builder_h(n, x);
                        vexpt.KO = vexpt.KO.setWp(Wp);

                        for i5 = m(m < i4)
%                              DispCounter(i5,i4-1,'i5:');

                            vexpt.nopre = vexpt.nopre.setFp(prm(i5),2);

                            comps(i1,i2,i2b,i3,i3b,i4,i5) = vexpt.InitRComp_top();
                            if comps(i1,i2,i3,i4,i5) > 0
                                y1 = BaselineWt(@CascadeBuilder, 8, prm(i1), [prm(i2b) prm(i2)], prm(i4));
                                y2 = BaselineWt(@CascadeBuilder, 8, x, [prm(i3b) prm(i3)], (prm(i4)));
                                wh = [wh; i1, i2, i2b, i3, i3b, i4, i5, x, y1 - y2, comps(i1,i2,i2b,i3,i3b,i4,i5)];
                            end
                        end%for i5
%                          DispCounter(i4,i4-1,'i5:');
                    end%for i4
                    DispCounter(M+1,M,'i4:');
                end%for i3b
                DispCounter(M+1,M,'i3b:');
            end%for i3
            DispCounter(M,M-1,'i3:');
        end%for i2b
        DispCounter(M,M-1,'i2b:');
    end%for i2
    DispCounter(M-1,M-2,'i2:');
end%for i1
DispCounter(M+1,M,'i1:');

end

