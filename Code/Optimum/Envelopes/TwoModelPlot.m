function TwoModelPlot( chains,twochains,ind,varargin )
%TWOMODELPLOT Summary of this function goes here
%   Detailed explanation goes here


srange=[twochains.s];
loglog(1./srange,srange.*[chains.A],'g',1./srange,[twochains.A].*srange,'r',1/twochains(1).sc,twochains(1).sc*twochains(1).Ac,'rd',varargin{:});
hold on
yl=ylim;
if exist('ind','var') && ~isempty(ind)
    loglog(1./srange,twochains(ind).snrb,'b',varargin{:});
    line([1 1]/srange(ind),yl,'Color','k',varargin{:});
end
hold off;

end

