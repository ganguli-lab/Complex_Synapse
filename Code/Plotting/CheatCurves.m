function [ s,env ] = CheatCurves( Wlen,t,x,varargin )
%S=CHEATCURVES(W,T,X) Encode with triangular decomposition of W, forget
%with shortcut removed/reduced version of W;
%S=CHEATCURVES(LEN,T,X) W randomly generated.
%   LEN = # states,
%   W = forgetting process
%   T = time values
%   x = fraction of shortcut reduction, 


persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='istransient';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addRequired('Wlen',@(x)validateattributes(x,{'numeric'},{'2d','square'},'CheatCurves','W/len',1))
    p.addRequired('t',@(x)validateattributes(x,{'numeric'},{'row'},'CheatCurves','t',2))
    p.addRequired('x',@(x)validateattributes(x,{'numeric'},{'row','nonnegative','<=',1},'CheatCurves','x',3))
    p.addOptional('sparsity',1,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1},'CheatCurves','sparsity',4));
    p.addParameter('cheatq',true,@(x) validateattributes(x,{'logical'},{'scalar'},'CheatCurves','cheatq'));
    p.addParameter('cheatW',true,@(x) validateattributes(x,{'logical'},{'scalar'},'CheatCurves','cheatW'));
end
p.parse(Wlen,t,x,varargin{:});
r=p.Results;
t=r.t;
x=r.x;


if isscalar(r.Wlen)
    len=r.Wlen;
    error(CheckValue(len,@(x) mod(x,2)==0,'even'));%even
    W=RandTrans(len,r.sparsity);
else
    W=r.Wlen;
    error(CheckSize(W,@(x) mod(length(x),2)==0,'even dim'));%even dim
    error(CheckValue(W,@isstochastic));
    len=length(W);
end


w=[-ones(len/2,1);ones(len/2,1)];
fp=0.5;

[Wp,Wm,w]=TriangleDcmp(W,fp,w);
q=Wp-Wm;
W=Wm+fp*q;

s=nan(length(x),length(t));


for i=1:length(x)
    if r.cheatq
        newq=(1-x(i))*q + x(i)*CheatInitWq(W,q,w);
    end
    if r.cheatW
        newW=RemoveShortcuts(W,x(i));
        if any(any(isnan(newW) | isinf(newW))) || abs(rcond(ones(size(newW))-newW))<1e-15
            continue
        end
    end
    if r.cheatq && r.cheatW
        s(i,:)=SNRcurveWq(t,newW,newq,fp,w);
    elseif r.cheatW
        s(i,:)=SNRcurveWq(t,newW,q,fp,w);
    elseif r.cheatq
        s(i,:)=SNRcurveWq(t,W,newq,fp,w);
    else
         s(i,:)=SNRcurveWq(t,W,q,fp,w);
   end%if cheat
end%for i

env=SNRenvelope(t,len);

if nargout==0
    plot(t,real([s;env]),varargin{:});
    legend([cellstr(num2str(x')); {'envelope'}]);
end%if nargout

% if any(diag(RemoveShortcuts(W))<-1)
%     disp('bad');
% end

end

