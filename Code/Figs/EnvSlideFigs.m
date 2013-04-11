fig=figure('PaperPositionMode','auto','Position',[60 60 800 400]);
ax=axes('Parent',fig,'FontSize',16);
hold(ax,'on')

t=0.5:0.5:4000;
n=40;

[env,bnds]=SNRenvelope(t,n);
[env2,bnds2]=SNRenvelope2(t,n);
yl=[env(end)/5 env2(1)];


t1=200;
s1=OneTimeMax2(t,t1,n);
tn1=find(t==t1);

lh1=line([1;1]*t1,yl','Color','r','LineWidth',2,'LineStyle','--','Parent',ax);
ph1=plot(t,s1,'b','LineWidth',2,'Parent',ax);
mh1=plot(t1,s1(tn1),'dg','Parent',ax,'MarkerFaceColor','g');



set(ax,'XScale','log','YScale','log')
xlim(ax,[t(1) t(end)])
ylim(ax,yl)
xlabel(ax,'Time','FontSize',30)
ylabel(ax,'SNR','FontSize',30);

pause;
print(fig,'env2_t1.eps','-depsc');

delete([lh1,ph1]);
clear('t1','s1','tn1','lh1','ph1');

t2=20;
s2=OneTimeMax2(t,t2,n);
tn2=find(t==t2);

lh2=line([1;1]*t2,yl','Color','r','LineWidth',2,'LineStyle','--','Parent',ax);
ph2=plot(t,s2,'b','LineWidth',2,'Parent',ax);
mh2=plot(t2,s2(tn2),'dg','Parent',ax,'MarkerFaceColor','g');

pause;
print(fig,'env2_t2.eps','-depsc');


delete([mh1,mh2,lh2,ph2]);
clear('mh1','mh2','t2','s2','tn2','lh2','ph2');

plot(t,env2,'g','LineWidth',3,'Parent',ax);
lh2=line([1;1]*bnds2',yl'*ones(1,length(bnds2)),'Color','k','LineWidth',1.5,'LineStyle','--','Parent',ax);

pause;
print(fig,'env2_all.eps','-depsc');

nvals=[4 10 16 28 40];
epsvals=0.1.^(1:0.5:3);
for nn=nvals
s=UniSMScurve(t,nn,1);
plot(t,s,'b','LineWidth',1.5)
end
for eps=epsvals
s=UniSMScurve(t,n,1,eps);
plot(t,s,'b','LineWidth',1.5)
end

pause;
print(fig,'env2_ex.eps','-depsc');

delete(lh2);
clear('nvals','epsvals','nn','eps','s','lh2');
plot(t,env2,'g',t,env,'g--','LineWidth',3,'Parent',ax);
line([1;1]*bnds',yl'*ones(1,length(bnds)),'Color','k','LineWidth',1.5,'LineStyle','--','Parent',ax);

pause;
print(fig,'env23.eps','-depsc');

close(fig);