function s=UniSMScurve(t,n,q,varargin)

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='UniSMScurve';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addRequired('t',@(x)validateattributes(x,{'numeric'},{'row'},'UniSMScurve','t',1))
    p.addRequired('n',@(x)validateattributes(x,{'numeric'},{'scalar','even'},'UniSMScurve','n',2))
    p.addRequired('q',@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1'},'UniSMScurve','q',3))
    p.addOptional('eps',1,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1'},'UniSMScurve','eps',4));
end
p.parse(t,n,q,varargin{:});
r=p.Results;
if any(strcmp('eps',p.UsingDefaults))
    r.eps=r.q;
end

qq=r.q*ones(1,r.n-1);
qq(1)=r.eps;
[Wp,Wm,w]=MakeSMS(qq);
s=SNRcurve(r.t,Wp,Wm,0.5,w);















end