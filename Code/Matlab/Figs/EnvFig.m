t=0.5:1:4000;
PlotEnvs( t,40 );
nvals=[4 10 16 28 40];
for n=nvals
s=DiffJumpcurve(t,n);
plot(t,s,'b','LineWidth',1.0)
end
epsvals=0.1.^(1:0.5:3);
for eps=epsvals
s=UniSMScurve(t,12,1,eps);
plot(t,s,'b','LineWidth',1.0)
end
legend({'Envelope','Conjectured envelope','Examples'},'Location','NorthEast')
set(gcf, 'PaperPositionMode', 'auto');
hold off
clear t s
clear n eps nvals epsvals