function [ h ] = PlotSMScurve( t,qv,fp,varargin )
%h=PLOTSMSCURVE(t,qv,fp,...) Plot SNR curve for Symmetric-Multistate model
%   t  = Time values
%   qv = Transition probabilities (row vector)
%   fp = Fraction of potentiation transitions (scalar in[0,1])
%   ...  Passed to 'plot'

error(CheckSize(qv,@isrow));
error(CheckSize(qv,@isrow));
error(CheckSize(fp,@isscalar));
error(CheckValue(fp,@(x) inrange(fp,0,1),'inrange(0,1)'));%fp in [0,1]


[Wp,Wm,w]=MakeSMS(qv);
s=SNRcurve(t,Wp,Wm,fp,w);
h=plot(t,s,varargin{:});

end

