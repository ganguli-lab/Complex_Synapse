function [ h ] = AreaCoeffPertPlot( Wp, Wm, Pertp, Pertm, fp, pertvals, varargin )
%H=AREACOEFFPERTPLOT(WP,WM,PERTP,PERTM,FP,PERTVALS,...) plot of area
%coefficients for a perturbation (see AreaCoeff)
%   H=plot handle
%   WP,WM = base transition rates for potentiation,depression
%   PERTP,PERTM = Perturbation of transition rates for potentiation,depression
%   FP = fraction of potentiation events
%   PERTVALS = vector/matrix of perturbation multipliers for plot
%   ... = passed to plot function

assert(ismat(Wp));%matrix
assert(issquare(Wp));%square
assert(samesize(Wp,Wm));%also square matrix of same size
assert(samesize(Wp,Pertp));%also square matrix of same size
assert(samesize(Wp,Pertm));%also square matrix of same size
assert(isscalar(fp));
assert(0<=fp && fp<=1);%fp in [0,1]
assert(isrow(pertvals));

% W=fp*Wp+(1-fp)*Wm;
% dW=fp*Pertp+(1-fp)*Pertm;

cc=zeros(size(Wp,1),length(pertvals));


for i=1:numel(pertvals)
%     cc(:,i)= AreaCoeffSorted(Wp + pertvals(i) * Pertp, Wm + pertvals(i) * Pertm, fp,w);
    cc(:,i)= AreaCoeff(Wp + pertvals(i) * Pertp, Wm + pertvals(i) * Pertm, fp);
end

h=plot(pertvals,cc,varargin{:});

end

