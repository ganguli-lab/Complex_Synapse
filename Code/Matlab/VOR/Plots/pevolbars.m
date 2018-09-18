function [ h ] = pevolbars( pevol,tind )
%h=PEVOLBARS(pevolbars) Plot time slice of evolution of state ditribution
%   For slides.
%   pevol = output of VORexptKO.ProbEvolsData (time x state)
%   tind  = which time slice?

yl=[0 0.2];
yt=yl;
h=bar(pevol(tind,:));
xlim([0.5 10.5])
ylim(yl);
set(gca,'YTick',yt);

print(gcf,['learn_anim_' int2str(tind) '.eps'],'-depsc')


end

