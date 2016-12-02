%%
builder_h = @SerialBuilder;
n = 10;
pwt = [0.12 0.14];
pko = 0.2;
df = -0.39;
df = {0.5, 0.5+df, 0.5-df};
tpre = 100;
train=150;
%%
builder_h = @CascadeBuilder;
n = 10;
pwt = [0.386 0.398];
pko = 0.466;
df = {0.522, 0.37, 0.998};
tpre = 200;
train=500;
%%
builder_h = @NonuniBuilder;
n = 10;
pwt = [0.4 0.4];
pko = 0.53;
df = {0.5, 0.3, 0.9};
tpre = 500;
train=800;
%%
vexpt=VORbuilderKO(builder_h, n, pwt(1), pwt(2), pko, df{:},train,tpre, false);
vexpt.PlotLearnS('all',false,'LineWidth',2);
xlim([0 train]);
%%
title('Serial model','FontSize',vexpt.FontSize);
print('serial_longterm.eps','-depsc');
%%
title('Cascade model','FontSize',vexpt.FontSize);
print('cascade_longterm.eps','-depsc');
%%
title('Non-uniform multistate model','FontSize',vexpt.FontSize);
print('nonuni_longterm.eps','-depsc');