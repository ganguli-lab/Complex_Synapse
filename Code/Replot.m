[Wp,Wm,w]=MakeSMS(q);
s=SNRcurve(t,Wp,Wm,0.5,w);
delete(h);
h=plot(t,s);