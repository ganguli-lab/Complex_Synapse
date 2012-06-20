t=0.1:0.1:4000;
n=20;
x=0.5;
beta=0.25;
%%
q=(x.^abs((1:(n-1))-n/2))/(1-x);
[Wp,Wm,w]=CascadeMSinterp(q,q,0.5,1);
casc=SNRcurve(t,Wp,Wm,0.5,w);
%%
[Wp,Wm,w]=MakeSMS((1-x)*q);
m1=SNRcurve(t,Wp,Wm,0.5,w);
%%
qm=ones(1,n-1);
qm(1:2:end)=beta;
[Wp,Wm,w]=MakeSMS(qm);
m2=SNRcurve(t,Wp,Wm,0.5,w);
%%
plot(t,casc,'b',t,m1,'g',t,m2,'r','LineWidth',2)
set(gca,'XScale','log','YScale','log')
xlim([t(1) t(end)])
ylim([m1(end) casc(1)])
xlabel('Time')
ylabel('SNR')
legend({'Cascade','Casc MS','Alt MS'},'Location','Best')
embiggen
set(gcf, 'PaperPositionMode', 'auto');
%%
clear Wp Wm w x beta n