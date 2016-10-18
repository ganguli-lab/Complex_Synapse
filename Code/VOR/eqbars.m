vexpt = VORbuilder(@SerialBuilder, 10, 0.12, 0.14, 0.2, 0.5, 0.11, 0.89, 5, 100, false);
yl=[0 0.2];
yt=yl;
xl=[0.5 10.5];

mem=vexpt.WT;
mem=mem.setFp(vexpt.withpre.fps(1));
p=mem.EqProb;
bar(p)
xlim(xl)
ylim(yl)
set(gca,'YTick',yt)
print(gcf,'serial_bar_wt_wo.eps','-depsc')

mem=mem.setFp(vexpt.withpre.fps(2));
p=mem.EqProb;
bar(p)
xlim(xl)
ylim(yl)
set(gca,'YTick',yt)
print(gcf,'serial_bar_wt_w.eps','-depsc')

mem=vexpt.KO;
mem=mem.setFp(vexpt.withpre.fps(1));
p=mem.EqProb;
bar(p)
xlim(xl)
ylim(yl)
set(gca,'YTick',yt)
print(gcf,'serial_bar_ko_wo.eps','-depsc')

mem=mem.setFp(vexpt.withpre.fps(2));
p=mem.EqProb;
bar(p)
xlim(xl)
ylim(yl)
set(gca,'YTick',yt)
print(gcf,'serial_bar_ko_w.eps','-depsc')

%%

vexpt = VORbuilder(@SerialBuilder, 2, 0.3, 0.3, 0.4, 0.5, 0.2, 0.8, 5, 100, false);
yl=[0 0.2];
yt=yl;
xl=[0.5 10.5];

mem=vexpt.WT;
mem=mem.setFp(vexpt.withpre.fps(1));
p=mem.EqProb;
bar(p)
xlim(xl)
ylim(yl)
set(gca,'YTick',yt)
print(gcf,'binary_bar_wt_wo.eps','-depsc')

mem=mem.setFp(vexpt.withpre.fps(2));
p=mem.EqProb;
bar(p)
xlim(xl)
ylim(yl)
set(gca,'YTick',yt)
print(gcf,'binary_bar_wt_w.eps','-depsc')

mem=vexpt.KO;
mem=mem.setFp(vexpt.withpre.fps(1));
p=mem.EqProb;
bar(p)
xlim(xl)
ylim(yl)
set(gca,'YTick',yt)
print(gcf,'binary_bar_ko_wo.eps','-depsc')

mem=mem.setFp(vexpt.withpre.fps(2));
p=mem.EqProb;
bar(p)
xlim(xl)
ylim(yl)
set(gca,'YTick',yt)
print(gcf,'binary_bar_ko_w.eps','-depsc')

