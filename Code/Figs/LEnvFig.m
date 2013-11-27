% nrange=[2 4 8 12];
nrange=2:2:12;
cla;
eh=loglog(1./srange,10*A.*srange,'g','LineWidth',3);
xlabel('Time');ylabel('Running average SNR');
embiggen;
xlim([1/max(srange) 1/min(srange)]);
set(gca,'XTick',[1e-3,1e-2,1e-1,1,1e1,1e2,1e3]);
hold on;
%
[Wp,Wm,w]=SerialBuilder(nrange(end),1);
%     ch=loglog(1./srange,10*SNRrunAve(1./srange,Wp,Wm,0.5,w),'r:','LineWidth',1.5);
uh=loglog(1./srange,10*SNRrunAve(1./srange,Wp,Wm,0.5,w),'b-','LineWidth',2);
%
for n=nrange(1:end-1)
%     [Wp,Wm,w]=DiffJump(n);
%     sh=loglog(1./srange,10*SNRrunAve(1./srange,Wp,Wm,0.5,w),'b--','LineWidth',1.5);
    [Wp,Wm,w]=SerialBuilder(n,1);
%     ch=loglog(1./srange,10*SNRrunAve(1./srange,Wp,Wm,0.5,w),'r:','LineWidth',1.5);
    ch=loglog(1./srange,10*SNRrunAve(1./srange,Wp,Wm,0.5,w),'b--','LineWidth',2);
end
epsvals=0.1.^(0.5:0.5:2);
for eps=epsvals
    q=ones(1,11);
    q(1)=eps;
    [Wp,Wm,w]=MakeSMS(q);
    sh=loglog(1./srange,10*SNRrunAve(1./srange,Wp,Wm,0.5,w),'b:','LineWidth',2);
end
% legend([eh ch sh],{'Numeric Envelope','Serial','Shortcut'},'Location','Best');
legend([eh uh ch sh],{'Numeric Envelope','Uniform serial','Shorter chains','Sticky end states'},'Location','Best');
clear nrange eh sh ch Wp Wm w n;
clear epsvals eps q uh;
xlim([1e-2 1e3])