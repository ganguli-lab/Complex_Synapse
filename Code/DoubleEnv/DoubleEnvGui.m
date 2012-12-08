function DoubleEnvGui( t1,t2,S2,n )
%DOUBLEENVGUI Summary of this function goes here
%   Detailed explanation goes here

[ h ] = PlotDoubleEnv(t1,t2,S2,n,[],[],[],[],[]);
axes(h(2));
set(h(2),'ButtonDownFcn',@replot);


    function replot(~,~)
        [t2,S2]=ginput(1);
        h=PlotDoubleEnv(t1,t2,S2,n,h(1),h(2),h(3),h(4),h(5));
    end



end

