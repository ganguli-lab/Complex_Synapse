function [ KTp,KTm ] = KTmults( Wp,Wm,dWp,dWm,varargin )
%[KTp,KTm]=KTMULTS(Wp,Wm,dWp,dWm,gamma) Kuhn-Tucker multipliers for maximisation
%of objective function wrt elements of Wp,Wm
%   Wp    = potentiation transition rates
%   Wm    = depression transition rates
%   dWpm  = gradient of objective function wrt elements of Wpm
%   gamma = size of diagonal constraint (default=1)
%   Assumes Wpm(ii) = -sum_(j~=i) Wpm(ij) enforced during all variations
%   Constraints are: W(ij)>=0 for i~=j
%                    W(ii)>=-gamma;

if nargin>4
    gamma=varargin{1};
else
    gamma=1;
end

mupdg = sum( Wp .* dWp, 2)/gamma;
mumdg = sum( Wm .* dWm, 2)/gamma;

ev=ones(size(mupdg))';

KTp = mupdg*ev - dWp;
KTm = mumdg*ev - dWm;

end

