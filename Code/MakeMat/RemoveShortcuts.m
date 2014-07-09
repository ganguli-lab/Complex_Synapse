function [ Wnsh ] = RemoveShortcuts( W, varargin )
%W_nsh=REMOVESHORTCUTS(W_sh,x) Remove shortcuts with perturbation that
%preserves equilibrium probabilities.
%   W_sh = transition rates with shortcuts
%   x = fraction of shortcut reduction, default=1
%   W_nsh = transition rates with no/reduced shortcuts

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='RemoveShortcuts';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addRequired('W',@(x)validateattributes(x,{'numeric'},{'2d','square'},'RemoveShortcuts','W',1))
    p.addOptional('x',1,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1},'RemoveShortcuts','x',2));
end
p.parse(W,varargin{:});
r=p.Results;

p=EqProb(r.W);
Wnsh=r.W;

for i=1:(length(r.W)-2)
    for j=(i+2):length(r.W)
        [dpertp,dpertm]=ShortcutPert(r.W,i,j-i);
        Wnsh=Wnsh - r.x * (p(i)*Wnsh(i,j)*dpertp + p(j)*Wnsh(j,i)*dpertm);
    end%for j
end%for i

end

