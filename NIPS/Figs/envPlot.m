h=plot(t,10*env,'g',tn,10*envnum,'r--',t,10*envshort,'r','LineWidth',3);
set(gca,'XScale','log','YScale','log');
ylim(yl);
xlim([t(1) t(end)]);
xlabel('Time');
ylabel('SNR');
line([1;1]*bnds',yl'*ones(1,length(bnds)),'Color','k','LineWidth',1.5,'LineStyle','--');
embiggen(gca,20);
legend({'envelope','numerical search','hand designed'},'Location','Best');

txFontSize=14;

    [x,y]=dsxy2figxy([t(1) bnds], (yl(1)^0.96*yl(2)^0.04)*[1 1]);
    annotation('doublearrow',x,y);
    pos=dsxy2figxy([t(1)^0.9*bnds^0.1 (yl(1)^0.95*yl(2)^0.05) t(1)^0.1*bnds^0.9 (yl(1)^0.9*yl(2)^0.1)]);
    annotation('textbox',pos,'String','Initial SNR bound active','VerticalAlignment','top','LineStyle','none','HorizontalAlignment','center','FontSize',txFontSize);
    
    [x,y]=dsxy2figxy([t(1) t(end)], (yl(1)^0.89*yl(2)^0.11)*[1 1]);
    annotation('doublearrow',x,y);
    pos=dsxy2figxy([t(1)^0.9*bnds^0.1 (yl(1)^0.88*yl(2)^0.12) t(1)^0.1*bnds^0.9 (yl(1)^0.83*yl(2)^0.17)]);
    annotation('textbox',pos,'String','Area bound active','VerticalAlignment','bottom','LineStyle','none','HorizontalAlignment','center','FontSize',txFontSize);
