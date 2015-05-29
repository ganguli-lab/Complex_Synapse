function [ fig ] = RunAveProven( srange,M,sqrtNum )
%RUNAVEPROVEN Summary of this function goes here
%   Detailed explanation goes here

fig=figure;
ax=axes('Parent',fig);

tau=wrev(1./srange);

AenvProven=(sqrtNum)*(M-1)./(tau+(M-1));

loglog(ax, tau,AenvProven,'g',...
    'LineWidth',2);
ylim([AenvProven(1)*srange(1) 5*(sqrtNum)]);
xlabel('Mean recall time, $\tau$','Interpreter','latex','FontSize',16);
ylabel('Recognition performance, $\overline{\mathrm{SNR}}(\tau)$','Interpreter','latex','FontSize',16);
title('Proven envelope','Interpreter','latex','FontSize',20);

area = (sqrtNum)*(M-1)./tau;
initial = (sqrtNum)*ones(size(tau));

toend=find(tau > (M-1)^0.9 * tau(end)^0.1, 1,'first');
fromstart=find(tau < (M-1)^0.5 * tau(1)^0.5, 1,'last');

hold(ax, 'on');
loglog(ax, tau(1:fromstart), initial(1:fromstart), 'k--',...
    tau(toend:end), area(toend:end), 'k--',...
    'LineWidth',1);


annotation('textbox',...
    dsxy2figxy(ax,[tau(1) sqrtNum tau(fromstart) 4*sqrtNum]),...
    'String','$\sqrt{N}$',...
    'Interpreter','latex','FontSize',16,...
    'LineStyle','none','VerticalAlignment','bottom','HorizontalAlignment','center');

annotation('textbox',...
    dsxy2figxy(ax,[tau(toend) 0.5*area(toend) tau(end) sqrt(area(toend)*5*sqrtNum)]),...
    'String','$\frac{\sqrt{N}(M-1)}{r\tau}$',...
    'Interpreter','latex','FontSize',22,...
    'LineStyle','none','VerticalAlignment','bottom','HorizontalAlignment','center');

annotation('textbox',...
    dsxy2figxy(ax,[sqrt((M-1)/srange(end)) 0.1*sqrtNum (M-1) 0.6*sqrtNum]),...
    'String','$\frac{\sqrt{N}(M-1)}{r\tau+(M-1)}$',...
    'Interpreter','latex','FontSize',22,...
    'LineStyle','none','VerticalAlignment','top','HorizontalAlignment','right');



end

