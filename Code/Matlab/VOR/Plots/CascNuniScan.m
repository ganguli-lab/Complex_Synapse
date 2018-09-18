%cascade/nonuni?
cn = true;
%%
cn = false;
%%
%parameter ranges
paramvals = (0.05:0.1:0.95);
ns = 4:10;
reps=20;
if cn
    ns = 4:2:16;
    paramvals = (0.05:0.05:0.5);
end
%%
%calculations
ncomps=ScanCNNum(paramvals,ns,reps,cn,'StepTolerance',1e-12,'ConstraintTolerance',1e-12);
%%
%plot
txopts={'Interpreter','latex','FontSize',20};
ax=gca;
bh=bar(ax, ns, ncomps','stacked');
% bh=bar(ax, ns, ncomps(:,1:6)','stacked');
ax.FontSize=16;
bh(1).FaceColor='none';
bh(1).EdgeColor='none';
bh(2).FaceColor=[0 0.4470 0.7410];
xlabel('$\#$ of synapstic states', txopts{:});
ylabel('max/min $\{\dot{L}_{\mathrm{WT}}(0) - \dot{L}_{DKO}(0)\}$', txopts{:});
if cn
    title('Cascade model, no pre', txopts{:});
    xlim([ns(1)-1 ns(end)+1]);
    ylim([-2 -1e-12]);
else
    title('Nonuniform multistate model, no pre', txopts{:});
    xlim([ns(1)-0.5 ns(end)+0.5]);
    ylim([-1 -1e-7]);
end
ax.YAxis.Scale = 'log';
%%
if cn
    print('cascade_scan.eps','-depsc');
else
    print('nonuni_scan.eps','-depsc');
end
