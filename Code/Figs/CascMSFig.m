t=10.^(-1:0.2:3);
x=0.5;
q=(x.^abs((1:19)-10))/(1-x);
[Wp,Wm,w]=CascadeMSinterp(q,q,0.5,1);
casc=SNRcurve(t,Wp,Wm,0.5,w);
%%
qm=ones(1,19);
[Wp,Wm,w]=MakeSMS(qm);
m1=SNRcurve(t,Wp,Wm,0.5,w);
%%
plot(t,casc,'b',t,m1,'g','LineWidth',2)
set(gca,'XScale','log','YScale','log')
xlim([t(1) t(end)])
ylim([casc(end) casc(1)])
xlabel('Time')
ylabel('SNR')
legend({'Cascade','Serial'},'Location','Best')
embiggen(gca,20)
set(gcf, 'PaperPositionMode', 'auto');
%%
clear x q Wp Wm w t m1 casc