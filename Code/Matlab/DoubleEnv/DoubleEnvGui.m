function DoubleEnvGui( t1,t2,S2,n )
%DOUBLEENVGUI Summary of this function goes here
%   Detailed explanation goes here

[ hs ] = PlotDoubleEnv(t1,t2,S2,n,[]);
axes(hs.ax_snr);
set(hs.ax_snr,'ButtonDownFcn',@replot);


    function replot(~,~)
        [t2,S2]=ginput(1);
        if S2 <= SNRenvelope2(t2,n)
            hs=PlotDoubleEnv(t1,t2,S2,n,hs);
        else
            warndlg('Click below envelope','Bad click');
        end
    end



end

