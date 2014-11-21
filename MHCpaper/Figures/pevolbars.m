function [ h ] = pevolbars( pevol,tind )
%PEVOLBARS Summary of this function goes here
%   Detailed explanation goes here

yl=[0 0.2];
yt=yl;
bar(pevol(tind,:));
xlim([0.5 10.5])
ylim(yl);
set(gca,'YTick',yt);

print(gcf,['learn_anim_' int2str(tind) '.eps'],'-depsc')


end

