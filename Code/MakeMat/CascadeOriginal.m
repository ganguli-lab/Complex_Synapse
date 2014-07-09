function [Wp,Wm,w]= CascadeOriginal( xp,xm,n,varargin )
%[Wp,Wm,w]=CASCADEORIGINAL(xp,xm,n) Makes Cascade model, as originally
%designed
%   XP,XM  = ratio of nearby transition rates
%   n  = number of states (even)
%   WP = potentiation transition rates
%   WM = depression transition rates
%   w  = Weights of states (+/-1)

% error(CheckSize(xp,@isscalar));
% error(CheckValue(xp,@(x) inrange(x,0,0.5)));
% error(CheckSize(xm,@isscalar));
% error(CheckValue(xm,@(x) inrange(x,0,0.5)));
% error(CheckSize(n,@isscalar));
% error(CheckValue(n,@isint));
persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='CascadeOriginal';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addRequired('xp',@(x)validateattributes(x,{'numeric'},{'scalar','positive'},'CascadeOriginal','xp',1))
    p.addRequired('xm',@(x)validateattributes(x,{'numeric'},{'scalar','positive'},'CascadeOriginal','xm',2))
    p.addRequired('n',@(x)validateattributes(x,{'numeric'},{'scalar','positive','integer'},'CascadeOriginal','n',3))
    p.addOptional('lambda',1,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1},'CascadeOriginal','lambda',4));
%     p.addParameter('Normalise',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
end
p.parse(xp,xm,n,varargin{:});
r=p.Results;

qp=(r.xp.^abs((1:r.n-1)-n/2))/(1-r.xp);
qm=(r.xm.^abs((1:r.n-1)-n/2))/(1-r.xm);
[Wp,~,w]=CascadeMSinterp(qp,qp,0.5,r.lambda);
[~,Wm]=CascadeMSinterp(qm,qm,0.5,r.lambda);

end

