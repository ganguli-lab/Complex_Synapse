% function PlotClique
% %PLOTCLIQUE Summary of this function goes here
% %   Detailed explanation goes here

n=zeros(1,100);
out=n;

for i=1:200
    n(i)=2*i+2;
    out(i)=max(CounterExample(n(i)));
    if mod(i,10)==0
        disp(int2str(i));
    end%if
end
plot(n,out);
% set(gca,'XScale','log','YScale','log')
% legend({'crq','init','cq','ev','Tmax','Eqbnd','L2','inittau','AI'});
hold off;
% end
clear i;
