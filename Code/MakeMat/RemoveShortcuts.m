function [ Wnsh ] = RemoveShortcuts( W, x )
%W_nsh=REMOVESHORTCUTS(W_sh,x) Remove shortcuts with perturbation that
%preserves equilibrium probabilities.
%   W_sh = transition rates with shortcuts
%   x = fraction of shortcut reduction, default=1
%   W_nsh = transition rates with no/reduced shortcuts

existsAndDefault('x',1);
error(CheckSize(W,@ismat));
error(CheckSize(x,@isscalar));
error(CheckValue(x,@(y) inrange(y,0,1),'inrange(0,1)'));

p=EqProb(W);
Wnsh=W;

for i=1:(length(W)-2)
    for j=(i+2):length(W)
        [dpertp,dpertm]=ShortcutPert(W,i,j-i);
        Wnsh=Wnsh - x * (p(i)*Wnsh(i,j)*dpertp + p(j)*Wnsh(j,i)*dpertm);
    end%for j
end%for i

end

