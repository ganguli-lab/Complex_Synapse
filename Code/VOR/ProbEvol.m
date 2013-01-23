function ProbEvol( Pt,t,titletext,varargin )
%PROBEVOL(Pt,t,titletext,...) Summary of this function goes here
%   Detailed explanation goes here

Parent=gca;
FontSize=10;
varargin=assignApplicable(varargin);
cla(Parent);

imagesc(t,1:size(Pt,2),Pt','Parent',Parent,varargin{:});
ylabel(Parent,'State');
xlabel(Parent,'Training time');
title(Parent,titletext);
h=colorbar('peer',Parent);
% set(get(h,'YLabel'),'String','Probability','Rotation',270,'VerticalAlignment','bottom','FontSize',FontSize);
colorbarlabel(h,'Probability','FontSize',FontSize);
        
end

