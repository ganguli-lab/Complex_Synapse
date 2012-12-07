function [ h ] = PlotDoubleEnv( t1,t2,S2,n )
%PLOTDOUBLEENV Summary of this function goes here
%   Detailed explanation goes here


fh=figure('Units','normalized','WindowStyle','docked');
ax_snr=axes('Parent',fh,'OuterPosition',[0,0.5,1,0.5]);
ax_exp=axes('Parent',fh,'OuterPosition',[0,0.32,1,0.16]);
ax_con1=axes('Parent',fh,'OuterPosition',[0,0.16,1,0.16]);
ax_con2=axes('Parent',fh,'OuterPosition',[0,0,1,0.16]);


S1=zeros(size(t1));
numexp=S1;
whichcase=S1;

for i=1:length(t1)
    [S1(i),whichcase(i),numexp(i)]=DoubleEnv(t1(i),t2,S2,n);
end

con1=mod(whichcase-1,2);
con2=mod((whichcase-1-con1)/2,2);

ph_snr=plot(t1,S1,'b',t2,S2,'r+','LineWidth',1.5,'Parent',ax_snr);
hold(ax_snr,'on')
[~,yl]=PlotEnvs(t1,n,'Parent',ax_snr,'format',false);
xlabel(ax_snr,'Time');
ylabel(ax_snr,'SNR');
ylim(ax_snr,yl);

ph_exp=area(t1,numexp,'Parent',ax_exp);
ylabel(ax_exp,'#exp');

ph_con1=area(t1,con1,'Parent',ax_con1);
ylabel(ax_con1,'SNR(0)');

ph_con2=area(t1,con2,'Parent',ax_con2);
ylabel(ax_con2,'Area');

h=[fh,ax_snr,ax_exp,ax_con1,ax_con2,ph_snr',ph_exp,ph_con1,ph_con2];

end

