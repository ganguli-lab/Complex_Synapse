t=10.^(-1:0.2:4);
x=0.5;
casc_model=SynapseMemoryModel.Build(@CascadeBuilder,0.5,{20,x});
casc=casc_model.SNRrunAve(t);
%%
serial_model=SynapseMemoryModel.Build(@SerialBuilder,0.5,{20,0.5});
serial=serial_model.SNRrunAve(t);
%%
plot(t,casc,'b',t,serial,'g','LineWidth',2)
set(gca,'XScale','log','YScale','log','FontSize',20)
xlim([t(1) t(end)])
ylim([casc(end) casc(1)])
% xlabel('Time','FontSize',20)
% ylabel('SNR','FontSize',20)
xlabel('Mean recall time, \tau','FontSize',20);
ylabel('$\overline{\mathsf{SNR}}\mathsf{(\tau)}$','Interpreter','latex','FontSize',20);
legend({'Cascade','Serial'},'Location','Best')
set(gcf, 'PaperPositionMode', 'auto');
%%
clear x q Wp Wm w t serial casc casc_model serial_model