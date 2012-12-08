function [ h ] = PlotDoubleEnv( t1,t2,S2,n,fh,ax_snr,ax_exp,ax_con1,ax_con2 )
%PLOTDOUBLEENV Summary of this function goes here
%   Detailed explanation goes here

if isempty(fh)
    fh=figure('Units','normalized','WindowStyle','docked');
end
if isempty(ax_snr)
    ax_snr=axes('Parent',fh,'OuterPosition',[0,0.5,1,0.5]);
end
if isempty(ax_exp)
    ax_exp=axes('Parent',fh,'OuterPosition',[0,0.32,1,0.16]);
end
if isempty(ax_con1)
    ax_con1=axes('Parent',fh,'OuterPosition',[0,0.16,1,0.16]);
end
if isempty(ax_con2)
    ax_con2=axes('Parent',fh,'OuterPosition',[0,0,1,0.16]);
end

cla(ax_snr);cla(ax_exp);cla(ax_con1);cla(ax_con2);

S1=zeros(size(t1));
numexp=S1;
whichcase=S1;

for i=1:length(t1)
    [S1(i),whichcase(i),numexp(i)]=DoubleEnv(t1(i),t2,S2,n);
end

con1=mod(whichcase-1,2);
con2=mod((whichcase-1-con1)/2,2);

xl=[t1(1) t1(end)];

[~,yl]=PlotEnvs(t1,n,'Parent',ax_snr,'format',false);
ph_snr=plot(t1,S1,'b',t2,S2,'r+','LineWidth',1.5,'Parent',ax_snr);
xlabel(ax_snr,'Time');
ylabel(ax_snr,'SNR');
ylim(ax_snr,yl);
xlim(ax_snr,xl);
set(ax_snr,'XScale','log','YScale','log');

ph_exp=area(t1,numexp,'Parent',ax_exp);
ylabel(ax_exp,'#exp');
xlim(ax_exp,xl);
set(ax_exp,'XScale','log');

ph_con1=area(t1,con1,'Parent',ax_con1);
ylabel(ax_con1,'SNR(0)');
xlim(ax_con1,xl);
set(ax_con1,'XScale','log');

ph_con2=area(t1,con2,'Parent',ax_con2);
ylabel(ax_con2,'Area');
xlim(ax_con2,xl);
set(ax_con2,'XScale','log');

h=[fh,ax_snr,ax_exp,ax_con1,ax_con2,ph_snr',ph_exp,ph_con1,ph_con2];

end

