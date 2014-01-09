yl=[0 0.2];
yt=yl;
mem=vexpt.WT;
mem=mem.setFp(vexpt.withpre.fps(1));
p=mem.EqProb;
bar(p)
xlim([0.5 10.5])
ylim(yl)
set(gca,'YTick',yt)
print(gcf,'serial_bar_wt_wo.eps','-depsc')
mem=mem.setFp(vexpt.withpre.fps(2));
p=mem.EqProb;
bar(p)
xlim([0.5 10.5])
ylim(yl)
set(gca,'YTick',yt)
print(gcf,'serial_bar_wt_w.eps','-depsc')
mem=vexpt.KO;
mem=mem.setFp(vexpt.withpre.fps(1));
p=mem.EqProb;
bar(p)
xlim([0.5 10.5])
ylim(yl)
set(gca,'YTick',yt)
print(gcf,'serial_bar_ko_wo.eps','-depsc')
mem=mem.setFp(vexpt.withpre.fps(2));
p=mem.EqProb;
bar(p)
xlim([0.5 10.5])
ylim(yl)
set(gca,'YTick',yt)
print(gcf,'serial_bar_ko_w.eps','-depsc')