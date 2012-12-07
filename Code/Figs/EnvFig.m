t=0.5:1:4000;
PlotEnvs( t,40 );
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
legend({'Envelope','Conjectured envelope','Examples'},'Location','NorthEast')
set(gcf, 'PaperPositionMode', 'auto');
hold off
clear t s