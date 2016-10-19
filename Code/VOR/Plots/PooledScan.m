%parameter ranges
paramvals = 0.05:0.1:0.95;
ns = 4:10;
%calculations
ncomps=ScanPooledNum(paramvals,4:10);
%%
%plot
txopts={'Interpreter','latex','FontSize',16};
ax=gca;
bh=bar(ax, ns-1, ncomps','stacked');
xlabel('Pool size', txopts{:});
ylabel('max/min learning rate difference', txopts{:});
title('Comparing WT: No pre $-$ w/ pre', txopts{:});
bh(1).FaceColor='none';
bh(1).EdgeColor='none';
bh(2).FaceColor=[0 0.4470 0.7410];
ax.YAxis.Scale = 'log';
%%
print('pooled_deponly_scan.eps','-depsc');
