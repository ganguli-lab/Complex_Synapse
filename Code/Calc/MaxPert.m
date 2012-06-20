function [ v ] = MaxPert( Wp, Wm, Pertp, Pertm, fp,w, pertvals )
%V=MAXPERT(WP,WM,PERTP,PERTM,FP,w,PERTVALS,...) value of perturbation that
%maximises area under SNR curve
%   V=value of maximising perturbation
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
assert(0<=fp && fp<=1);%fp in [0,1]

area=zeros(size(pertvals));

for i=1:numel(pertvals)
    area(i)= SNRarea(Wp + pertvals(i) * Pertp, Wm + pertvals(i) * Pertm, fp, w);
end

[~,I]=max(area);
v=pertvals(I);

end

