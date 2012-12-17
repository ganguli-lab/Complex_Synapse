function [ hs ] = PlotDoubleEnv( t1,t2,S2,n,hs,varargin )
%hs=PLOTDOUBLEENV(t1,t2,S2,n,hs) Plot Two-time envelope
%   t1 = time of maximisation
%   t2 = time of fixed SNR
%   S2 = SNR(t2)
%   n = # states
%   hs = handle struct, fields:
%       fh      = figure
%       ax_???  = axes for plot
%       ph_???  = plot handle
%       ??_snr  = envelopes
%       ??_exp  = # of exponentials for optimal curve
%       ??_conX = is constraint X active? (X=1,2)
%
%   Constraints:
%       1: sum_a c_a q_a < 1
%       2: sum_a c_a     < n-1
%       3: c_a sqrt(q_a) < gamma

Constraint3=false;
varargin=assignApplicable(varargin);

if isempty(hs)
    hs.fh=figure('Units','normalized','WindowStyle','docked');
    hs.ax_snr=axes('Parent',hs.fh,'OuterPosition',[0,0.5,1,0.5]);
    hs.ax_exp=axes('Parent',hs.fh,'OuterPosition',AxPos(1));
    hs.ax_con1=axes('Parent',hs.fh,'OuterPosition',AxPos(2));
    hs.ax_con2=axes('Parent',hs.fh,'OuterPosition',AxPos(3));
    if Constraint3
        hs.ax_con3=axes('Parent',hs.fh,'OuterPosition',Axpos(4));
    end
end

cla(hs.ax_snr);cla(hs.ax_exp);cla(hs.ax_con1);cla(hs.ax_con2);

S1=zeros(size(t1));
numexp=S1;
cons=zeros(2+Constraint3,length(t1));

for i=1:length(t1)
    [S1(i),whichcase,numexp(i)]=DoubleEnv(t1(i),t2,S2,n,'Constraint3',Constraint3,varargin{:});
    cons(:,i)=BinVec(whichcase,2,2+Constraint3)';
end



xl=[t1(1) t1(end)];

[hs.envh,yl]=PlotEnvs(t1,n,'Parent',hs.ax_snr,'Format',false);
hs.ph_snr=plot(t1,S1,'b',t2,S2,'r+','LineWidth',1.5,'Parent',hs.ax_snr);
xlabel(hs.ax_snr,'Time');
ylabel(hs.ax_snr,'SNR');
ylim(hs.ax_snr,yl);
xlim(hs.ax_snr,xl);
set(hs.ax_snr,'XScale','log','YScale','log');

hs.ph_exp=area(t1,numexp,'Parent',hs.ax_exp);
ylabel(hs.ax_exp,'#exp');
xlim(hs.ax_exp,xl);
set(hs.ax_exp,'XScale','log');

hs.ph_con1=area(t1,cons(1,:),'Parent',hs.ax_con1);
ylabel(hs.ax_con1,'Initial');
xlim(hs.ax_con1,xl);
set(hs.ax_con1,'XScale','log');

hs.ph_con2=area(t1,cons(2,:),'Parent',hs.ax_con2);
ylabel(hs.ax_con2,'Area');
xlim(hs.ax_con2,xl);
set(hs.ax_con2,'XScale','log');

if Constraint3
    hs.ph_con3=area(t1,cons(3,:),'Parent',hs.ax_con2);
    ylabel(hs.ax_con3,'Mode');
    xlim(hs.ax_con3,xl);
    set(hs.ax_con3,'XScale','log');
end

    function axpos=AxPos(n)
        numax=3+Constraint3;
        axpos=[0,0.5*(numax-n)/numax,1,0.5/numax];
    end


end

