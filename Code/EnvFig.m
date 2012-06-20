t=0.5:1:4000;
[env,bnds]=SNRenvelope(t,40);
[env2,bnds2]=SNRenvelope2(t,40);
plot(t,env2,'g',t,env,'g--','LineWidth',3);
hold on
set(gca,'XScale','log','YScale','log')
xlim([t(1) t(end)])
yl=[env(end) 1];
ylim(yl)
embiggen
xlabel('Time')
ylabel('SNR');
nvals=[4 10 16 28 40];
epsvals=0.1.^(1:0.5:3);
for n=nvals
s=UniSMScurve(t,n,1);
plot(t,s,'b','LineWidth',1.5)
end
for eps=epsvals
s=UniSMScurve(t,40,1,eps);
plot(t,s,'b','LineWidth',1.5)
end
line([1;1]*bnds',yl'*ones(1,length(bnds)),'Color','k','LineWidth',1.5,'LineStyle',':');
line([1;1]*bnds2',yl'*ones(1,length(bnds2)),'Color','k','LineWidth',1.5,'LineStyle','--');
legend({'Envelope','Conjectured envelope','Examples'},'Location','NorthEast')
set(gcf, 'PaperPositionMode', 'auto');
hold off
clear t env env2 bnds bnds2 yl s