% function PlotClique
% %PLOTCLIQUE Summary of this function goes here
% %   Detailed explanation goes here

n=zeros(1,100);
cq=n;
crq=n;
init=n;
ev=n;
Tmax=n;
Eqbnd=n;
L2=n;
AI=n;

for i=1:100
    n(i)=6*i;
    [ crq(i),init(i),cq(i),ev(i),Tmax(i),Eqbnd(i),L2(i),AI(i) ]=TestClique(n(i));
    if mod(i,10)==0
        disp(int2str(i));
    end%if
end
plot(n,real([crq;init;cq;ev;Tmax;Eqbnd;L2;init./sqrt(ev);AI]));
set(gca,'XScale','log','YScale','log')
legend({'crq','init','cq','ev','Tmax','Eqbnd','L2','inittau','AI'});
hold off;
% end
clear i;
