function ProbEvol( Pt,t,titletext,varargin )
%PROBEVOL(Pt,t,titletext,...) Summary of this function goes here
%   Detailed explanation goes here

Parent=gca;
FontSize=10;
axFontSize=10;
Pooled=false;
varargin=assignApplicable(varargin);
cla(Parent);

n=1:size(Pt,2);
if Pooled
    n=n-1;
    xlab='Number potentiated';
else
    xlab='State';
end


imagesc(t,n,Pt','Parent',Parent,varargin{:});
ylabel(Parent,xlab,'FontSize',axFontSize);
xlabel(Parent,'Training time','FontSize',axFontSize);
title(Parent,titletext,'FontSize',FontSize);
h=colorbar('peer',Parent);
% set(get(h,'YLabel'),'String','Probability','Rotation',270,'VerticalAlignment','bottom','FontSize',FontSize);
colorbarlabel(h,'Probability','FontSize',axFontSize);
        
end

