h=plot(t,10*env,'g',tn,10*envnum,'r--',t,10*envshort,'r','LineWidth',3);
set(gca,'XScale','log','YScale','log');
ylim(yl);
xlim([t(1) t(end)]);
xlabel('Time');
ylabel('SNR');
line([1;1]*bnds',yl'*ones(1,length(bnds)),'Color','k','LineWidth',1.5,'LineStyle','--');
legend({'envelope','numerical search','hand designed'},'Location','Best');