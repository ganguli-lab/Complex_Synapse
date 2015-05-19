function [ fig ] = RunAveProven( srange,M,sqrtNum )
%RUNAVEPROVEN Summary of this function goes here
%   Detailed explanation goes here

fig=figure;
ax=axes('Parent',fig);

AenvProven=(sqrtNum)*(M-1)./(1+(M-1)*srange);

loglog(ax, 1./srange,AenvProven.*srange,'g',...
    'LineWidth',2);
ylim([AenvProven(1)*srange(1) 5*(sqrtNum)]);
xlabel('Mean recall time, $\tau$','Interpreter','latex','FontSize',16);
ylabel('Recognition performance, $\overline{\mathrm{SNR}}(\tau)$','Interpreter','latex','FontSize',16);
title('Proven envelope','Interpreter','latex','FontSize',20);
area = (sqrtNum)*(M-1)*srange;
initial = (sqrtNum)*ones(size(srange));

toend=round(0.45*length(srange));
fromstart=round(0.4*length(srange));

hold(ax, 'on');
loglog(ax, 1./srange(toend:end), initial(toend:end), 'k--',...
    1./srange(1:fromstart), area(1:fromstart), 'k--',...
    'LineWidth',1);


annotation('textbox',...
    dsxy2figxy(ax,[1/srange(end) sqrtNum 1/srange(toend)-1/srange(end) 4*sqrtNum]),...
    'String','$\sqrt{N}$',...
    'Interpreter','latex','FontSize',16,...
    'LineStyle','none','VerticalAlignment','bottom','HorizontalAlignment','center');
annotation('textbox',...
    dsxy2figxy(ax,[1/srange(fromstart) sqrtNum 1/srange(1)-1/srange(fromstart) 4*sqrtNum]),...
    'String','$\frac{\sqrt{N}(M-1)}{r\tau}$',...
    'Interpreter','latex','FontSize',22,...
    'LineStyle','none','VerticalAlignment','bottom','HorizontalAlignment','center');



end

