function [ lengthmet,lengths,ix ] = ChainLength( qv,varargin )
%[lengthmet,lengths,ix]=CHAINLENGTH(qv,dim) estimate effective length of chain
%   qv  = nearest neighbour transitions
%   dim = dimension of qv corresponding to different states
%   lengthmet = metric measuring how much like a chain of that length it is
%   lengths   = length of chain corresponding to each row of lengthmet
%   ix        = where is lengthmet maximised (linear ind for row of lengthmet

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='ChainLength';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addRequired('qv',@(x)validateattributes(x,{'numeric'},{'2d'},'ChainLength','qv',1))
    p.addOptional('dim',1+isrow(qv),@(x)validateattributes(x,{'numeric'},{'scalar','positive','integer'},'ChainLength','dim',2));
end
p.parse(W,varargin{:});
qv=p.Results.qv;
dim=p.Results.dim;
if any(strcmp('dim',p.UsingDefaults))
    dim=1+isrow(qv);
end

n=size(qv,dim)+1;
qv=permute(qv,[dim 1:(dim-1) (dim+1):ndims(qv)]);

lengths=2:2:n;

siz=size(qv);
siz(1)=length(lengths);
lengthmet=zeros(siz(1),prod(siz(2:end)));

for i=1:siz(1)
    lengthmet(i,:)=min(qv(n/2:((n+lengths(i)-2)/2),:),[],1);
end

lengthmet=[-diff(lengthmet,1,1);lengthmet(end,:)];

[~,ix]=max(lengthmet,[],2);
ix(end)=find(lengthmet(end,:)>0.99,1,'last');

lengthmet=reshape(lengthmet,siz);


end

