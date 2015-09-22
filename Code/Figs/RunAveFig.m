
% Interpreter='latex';
Interpreter='tex';

load('laplacebnd12.mat','AenvAll','srange');
load('laplacebndChainA12.mat','AenvChains');

AenvProven=11./(1+11*srange);


loglog( 1./srange,100*AenvProven.*srange,'g',...
    1./srange,100*AenvAll.*srange,'r:',...
    1./srange,100*AenvChains.*srange,'r-',...
    'LineWidth',2);
%     1./srange,100*AenvHomC.*srange,'r--',...
ylim([1e-1 500]);
switch Interpreter
    case 'latex'
        xlabel('Mean recall time, $\tau$','Interpreter',Interpreter,'FontSize',16);
        ylabel('Recognition performance, $\overline{\mathrm{SNR}}(\tau)$','Interpreter',Interpreter,'FontSize',16);
    case 'tex'
        xlabel('Mean recall time, \tau','Interpreter',Interpreter,'FontSize',16);
        ylabel('\textsf{Recognition performance,} $\overline{\mathsf{SNR}}\mathsf{(\tau)}$','Interpreter','latex','FontSize',16);
end
title('Numerical envelopes','Interpreter',Interpreter,'FontSize',20);
legend({'Proven envelope',...
    'Numeric: all topologies',...
    'Numeric: serial topology',...
    },...
    'Location','southwest','Interpreter',Interpreter,'FontSize',16);