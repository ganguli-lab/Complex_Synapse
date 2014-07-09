function [ Wp,Wm,w ] = OneShortcut( n,i,j,varargin )
%[Wp,Wm,w]=ONESHORTCUT(n,i,j,qs,qd) Synapse with one shortcut for
%potentiation, flip symmetric for depression
%   n  = # states
%   i  = starting state for shortcut
%   j  = final state for shortcut
%   qs = transition prob of shortcut (default: 1)
%   qd = transition prob of direct link out of state i (default: 1-qs)

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='OneShortcut';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addRequired('n',@(x)validateattributes(x,{'numeric'},{'scalar','integer','positive'},'OneShortcut','n',1))
    p.addRequired('i',@(x)validateattributes(x,{'numeric'},{'scalar','integer','positive','<=',n},'OneShortcut','i',2))
    p.addRequired('j',@(x)validateattributes(x,{'numeric'},{'scalar','integer','positive','<=',n},'OneShortcut','j',3))
    p.addOptional('qs',1,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1},'OneShortcut','qs',4));
    p.addOptional('qd',0,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1},'OneShortcut','qd',5));
end
p.parse(n,i,j,varargin{:});
r=p.Results;
if any(strcmp('qd',p.UsingDefaults))
    r.qd=1-r.qs;
end

q=ones(1,r.n-1);
q(r.i)=r.qd;

[Wp,Wm,w]=MakeSMS(q);

Wp(r.i,r.j)=r.qs;
Wm(r.n-r.i+1,r.n-r.j+1)=r.qs;
Wp=StochastifyC(Wp);
Wm=StochastifyC(Wm);

end

