t=0:0.1:3;
s=5*exp(-t);
ts=log(5);
plot(t,s,'b','LineWidth',2)
line([0 ts ts],[1 1 0],'Color','g','LineStyle','--','LineWidth',3)
set(gca,'XTick',ts,'XTickLabel','lifetime','YTick',1,'YTickLabel','1')
xlabel('Time')
ylabel('SNR')
embiggen
set(gcf, 'PaperPositionMode', 'auto');
clear t s ts