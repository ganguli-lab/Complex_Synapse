function [ h,c ] = AreaCoeffBar( Wp, Wm, fp, w, varargin )
%H=AREACOEFFBAR(WP,WM,FP,w,...) Bar chart of area coeffs
%   WP = potentiation transition rates
%   WM = depression transition rates
%   FP = Fraction of potentiation transitions
%   w = Weights of states (+/-1)
%   H = handle
%   ... passed to bar

assert(ismat(Wp));%matrix
assert(issquare(Wp));%square
assert(samesize(Wp,Wm));%also square matrix of same size
assert(isscalar(fp));
assert(0<=fp && fp<=1);%fp in [0,1]
assert(iscol(w));%row
assert(length(w)==length(Wp));%same size
assert(all(abs(w)==1));%+/-1


c=AreaCoeff(Wp,Wm,fp);


cc=[c;c];
cc(1,:)=cc(1,:).*(1+w)'/2;
cc(2,:)=cc(2,:).*(1-w)'/2;


h=bar(cc',varargin{:});

Parent=gca;
varargin=assignApplicable(varargin);


legend(Parent,'w>0','w<0','Location','Best');

end

