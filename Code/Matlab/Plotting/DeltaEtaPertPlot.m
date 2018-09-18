function [ h ] = DeltaEtaPertPlot( Wp, Wm, Pertp, Pertm, fp, w, pertvals, varargin )
%H=DELTAETAPERTPLOT(WP,WM,PERTP,PERTM,FP,w,PERTVALS,...) plot of partial
%mixing times for a perturbation
%   H=plot handle
%   WP,WM = base transition rates for potentiation,depression
%   PERTP,PERTM = Perturbation of transition rates for potentiation,depression
%   FP = fraction of potentiation events
%   PERTVALS = vector/matrix of perturbation multipliers for plot
%   w  = Weights of states (+/-1)
%   ... = passed to plot function

assert(ismat(Wp));%matrix
assert(issquare(Wp));%square
assert(samesize(Wp,Wm));%also square matrix of same size
assert(samesize(Wp,Pertp));%also square matrix of same size
assert(samesize(Wp,Pertm));%also square matrix of same size
assert(isscalar(fp));
assert(0<=fp && fp<=1);%fp in [0,1]
assert(isrow(pertvals));
assert(iscol(w));%row
assert(length(w)==length(Wp));%same size
assert(all(abs(w)==1));%+/-1

W=fp*Wp+(1-fp)*Wm;
dW=fp*Pertp+(1-fp)*Pertm;

deta=zeros(size(Wp,1),length(pertvals));


for i=1:numel(pertvals)
    deta(:,i)= DeltaEta(W + pertvals(i) * dW,w);
end

h=plot(pertvals,deta,varargin{:});

end

