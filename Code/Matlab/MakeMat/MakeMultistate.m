function [ Wp,Wm,w ] = MakeMultistate( qp,qm )
%[WP,WM,w]=MAKEMULTISTATE(QP,QM) Transition rates for Multistate toplogy
%   QP/QM=adjacent transiton rates for potentiation/depression, vector

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='MakeMultistate';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addRequired('qp',@(x)validateattributes(x,{'numeric'},{'vector','nonnegative'},'MakeMultistate','qp',1))
    p.addRequired('qm',@(x)validateattributes(x,{'numeric'},{'vector','nonnegative'},'MakeMultistate','qm',2))
end
p.parse(qp,qm);
validateattributes(qm,{'numeric'},{'size',size(qp)},'MakeMultistate','qm',2)
r=p.Results;
error(CheckSize(r.qp,@(x) mod(length(x),2)==1,'odd dim'));

Wp=StochastifyC(diag(r.qp,1));
Wm=StochastifyC(diag(r.qm,-1));

if nargout>2
w=[-ones((length(qp)+1)/2,1);ones((length(qp)+1)/2,1)];
end

end

