%parameter ranges
paramvals = (0.05:0.1:0.95);
ns = 4:10;
reps=2;
%%
%calculations
ncomps=ScanPooledNum(paramvals,ns,'Top',reps);
%%
%plot
txopts={'Interpreter','latex','FontSize',20};
ax=gca;
bh=bar(ax, ns-1, ncomps','stacked');
xlabel('$\#$ of synapses in pool', txopts{:});
ylabel('max/min $\{\dot{L}_{\mathrm{WT}}(0) - \dot{L}_{\mathrm{K}^b\mathrm{D}^{b-/-}}(0)\}$', txopts{:});
title('Pooled resource model, WT', txopts{:});
bh(1).FaceColor='none';
bh(1).EdgeColor='none';
bh(2).FaceColor=[0 0.4470 0.7410];
ax.YAxis.Scale = 'log';
%%
print('pooled_scan.eps','-depsc');
