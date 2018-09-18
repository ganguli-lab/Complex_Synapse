function [ h ] = ScaledDerivPertPlot( Wp, Wm, Pertp, Pertm, fp, w, pertvals, varargin )
%H=AREAPERTPLOT(WP,WM,PERTP,PERTM,FP,w,PERTVALS,...) plot of area under SNR
%curve for a perturbation
%   H=plot handle
%   WP,WM = base transition rates for potentiation,depression
%   PERTP,PERTM = Perturbation of transition rates for potentiation,depression
%   FP = fraction of potentiation events
%   w = Weights of states (+/-1)
%   PERTVALS = vector/matrix of perturbation multipliers for plot
%   ... = passed to plot function

assert(ismat(Wp));%matrix
assert(issquare(Wp));%square
assert(samesize(Wp,Wm));%also square matrix of same size
assert(samesize(Wp,Pertp));%also square matrix of same size
assert(samesize(Wp,Pertm));%also square matrix of same size
assert(isscalar(fp));

darea=zeros(size(pertvals));

[row,col]=ind2sub(size(Pertp),find(Pertp==1));
% [row,col]=ind2sub(size(Pertm),find(Pertm==1));

for i=1:numel(pertvals)
    darea(i)= ScaledDeriv(Wp + pertvals(i) * Pertp, Wm + pertvals(i) * Pertm, fp, w, row,col);
end

h=plot(pertvals,darea,varargin{:});

end

