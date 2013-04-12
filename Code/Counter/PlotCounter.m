% function PlotClique
% %PLOTCLIQUE Summary of this function goes here
% %   Detailed explanation goes here
fig=figure('PaperPositionMode','auto','Position',[60 60 800 400]);
ax=axes('Parent',fig,'FontSize',16);
%allvars={'crqm','crq','c2q','init','cqm','ev','L2','AI','Eqbnd','Tmax'};
% varnames={'crq','crqm','c2q'};
%varnames={'crq','qPi','init','A','etaPi'};
varnames={'crqm','crq','c2q','AI','L2'};
numpts=100;
n=zeros(1,numpts);
output=zeros(length(varnames),numpts);

for i=1:numpts
    n(i)=2*i+4;
%     [ crq(i),init(i),cq(i),ev(i),Tmax(i),Eqbnd(i),L2(i),AI(i) ]=TestClique(n(i));
    output(:,i)=TestCounterEx(n(i),varnames);
    if mod(i,10)==0
        disp(int2str(i));
    end%if
end
plot(n,real(output),'LineWidth',3);
set(gca,'XScale','log','YScale','log')
varnames={'$\max_a I_a \sqrt{\tau_a}$','$\sum_a I_a \sqrt{\tau_a}$','$\sum_a I_a^2 \tau_a$','Area $\times$ SNR(0)','$\int \mathrm{SNR}(t)^2 \mathrm{d}t$'};
legend(varnames,'Interpreter','LaTeX','Location','Best');
xlabel('Number of states','FontSize',24);
xlim([n(1) n(end)]);
hold off;
% end
clear i numpts;
