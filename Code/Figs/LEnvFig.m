nrange=[2 4 8 12];
cla;
eh=loglog(1./srange,10*A.*srange,'g','LineWidth',3);
xlabel('Time');ylabel('Running average SNR');
embiggen;
xlim([1/max(srange) 1/min(srange)]);
set(gca,'XTick',[1e-3,1e-2,1e-1,1,1e1,1e2,1e3]);
hold on;
for n=nrange
    [Wp,Wm,w]=DiffJump(n);
    sh=loglog(1./srange,10*SNRrunAve(1./srange,Wp,Wm,0.5,w),'b--','LineWidth',1.5);
    [Wp,Wm,w]=SerialBuilder(n,1);
    ch=loglog(1./srange,10*SNRrunAve(1./srange,Wp,Wm,0.5,w),'r:','LineWidth',1.5);
end
legend([eh ch sh],{'Numeric Envelope','Serial','Shortcut'},'Location','Best');