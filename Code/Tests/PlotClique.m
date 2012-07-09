% function PlotClique
% %PLOTCLIQUE Summary of this function goes here
% %   Detailed explanation goes here

%allvars={'crqm','crq','c2q','init','cqm','ev','L2','AI','Eqbnd','Tmax'};
% varnames={'crq','crqm','c2q'};
varnames={'crq','Eqbnd','init','ev'};
numpts=100;
n=zeros(1,numpts);
output=zeros(length(varnames),numpts);

for i=1:numpts
    n(i)=2*i+4;
%     [ crq(i),init(i),cq(i),ev(i),Tmax(i),Eqbnd(i),L2(i),AI(i) ]=TestClique(n(i));
    output(:,i)=TestClique(n(i),varnames);
    if mod(i,10)==0
        disp(int2str(i));
    end%if
end
plot(n,real(output));
set(gca,'XScale','log','YScale','log')
legend(varnames);
hold off;
% end
clear i numpts;
