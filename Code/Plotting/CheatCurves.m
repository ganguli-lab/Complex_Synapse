function [ s,env ] = CheatCurves( Wlen,t,x,sparsity,varargin )
%S=CHEATCURVES(W,T,X) Encode with triangular decomposition of W, forget
%with shortcut removed/reduced version of W;
%S=CHEATCURVES(LEN,T,X) W randomly generated.
%   LEN = # states,
%   W = forgetting process
%   T = time values
%   x = fraction of shortcut reduction, 

error(CheckSize(t,@isrow));
error(CheckSize(x,@isrow));

cheatq=true;
cheatW=true;

varargin=assignApplicable(varargin);

existsAndDefault('sparsity',1);



if isscalar(Wlen)
    len=Wlen;
    error(CheckValue(len,@(x) mod(x,2)==0,'even'));%even
    W=RandTrans(len,sparsity);
else
    W=Wlen;
    error(CheckSize(W,@ismat));%matrix
    error(CheckSize(W,@issquare));%square
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
    if cheatq
        newq=(1-x(i))*q + x(i)*CheatInitWq(W,q,w);
    end
    if cheatW
        newW=RemoveShortcuts(W,x(i));
        if any(any(isnan(newW) | isinf(newW))) || abs(rcond(ones(size(newW))-newW))<1e-15
            continue
        end
    end
    if cheatq && cheatW
        s(i,:)=SNRcurveWq(t,newW,newq,fp,w);
    elseif cheatW
        s(i,:)=SNRcurveWq(t,newW,q,fp,w);
    elseif cheatq
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

