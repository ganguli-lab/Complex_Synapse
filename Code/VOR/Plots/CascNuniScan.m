%cascade/nonuni?
cn = true;
%%
%parameter ranges
paramvals = (0.05:0.1:0.95);
ns = 4:10;
reps=10;
if cn
    ns = 4:2:16;
    paramvals = paramvals/2;
end
%%
%calculations
ncomps=ScanCNNum(paramvals,ns,reps,cn,'FunctionTolerance',1e-7,'StepTolerance',1e-11,'OptimalityTolerance',1e-7);
%%
%plot
txopts={'Interpreter','latex','FontSize',20};
ax=gca;
bh=bar(ax, ns, ncomps','stacked');
xlabel('$\#$ of synapstic states', txopts{:});
ylabel('max/min $\{\dot{L}_{\mathrm{WT}}(0) - \dot{L}_{\mathrm{K}^b\mathrm{D}^{b-/-}}(0)\}$', txopts{:});
if cn
    title('Cascade model, no pre', txopts{:});
else
    title('Nonuniform multistate model, no pre', txopts{:});
end
bh(1).FaceColor='none';
bh(1).EdgeColor='none';
bh(2).FaceColor=[0 0.4470 0.7410];
% ax.YAxis.Scale = 'log';
%%
if cn
    print('cascade_scan.eps','-depsc');
else
    print('nonuni_scan.eps','-depsc');
end