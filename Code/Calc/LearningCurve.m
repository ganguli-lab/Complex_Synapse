function [ S ] = LearningCurve( W0,w,t,W1,varargin )
%S=LEARNINGCURVE(W0,w,t,W1,...) mean synaptic weight as function of time,
%starting in equilibrium state for W0, then evolving according to W1
%   S(t)=p(t)w
%   dp/dt = pW
%   p(0)W_0=0
%Additional arguments:
%   ...t_n,W_n+1,...
%       W=W_n for t_n-1 < t <= t_n
%       t_0=0,
%   actually uses max(t<t_n) in place of t_n.

error(CheckSize(W0,@isstochasticC));
error(CheckSize(w,@iscolumn));
error(CheckSize(t,@isrow));
error(CheckValue(t,@(x) all(diff(x)>0),'increasing'));
error(CheckSize(W1,@isstochasticC));

if mod(nargin,2)~=0
    error('Wrong number of arguments: each t_n must be followed by W_n+1');
end

numchange=(nargin-4)/2;

%p0: row=which state.
p0=EqProb(W0);
%V: col=which eigenmode, row=which state.
%D: col=row=which eigenmode
[V,D]=eig(W1);
%expqt,pt: col=which eigenmode, row=what time.
expqt=exp(diag(D)*t)';
pt=expqt*diag(p0*V)/V;
%S: col=what time.
S=(pt*w)';



for i=1:numchange
    tchange=varargin{2*i-1};
    W=varargin{2*i};
    
    error(CheckSize(tchange,@isscalar));
    error(CheckSize(W1,@isstochasticC));
    
    valid=t>tchange;
    ix=find(~valid,1,'last');
    
    p0=pt(ix,:);
    [V,D]=eig(W);
    expqt=exp(diag(D)*(t-t(ix)))';%row=wich eigenmode, col=what time.
    pt=(expqt*diag(p0*V))/V;%col=wich eigenmode, row=what time.
    
    newS=(pt*w)';
    S(valid)=newS(valid);
end


end

