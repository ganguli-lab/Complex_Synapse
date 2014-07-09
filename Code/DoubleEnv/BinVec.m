function [ binvec ] = BinVec( num,varargin )
%bv=BINREP(num,base,digits) Vector of digits in Binary represantation
%   num = number to represent
%   base = base of representation (default=2)
%   digits = number of digits to use
%num =sum_{i=1:digits} bv(i)*base^(i-1)

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='BinVec';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addRequired('num',@(x)validateattributes(x,{'numeric'},{'scalar'},'BinVec','W',1))
    p.addOptional('base',2,@(x)validateattributes(x,{'numeric'},{'scalar','positive'},'BinVec','base',2));
    p.addOptional('digits',0,@(x)validateattributes(x,{'numeric'},{'scalar','positive','integer'},'BinVec','digits',2));
end
p.parse(num,varargin{:});
r=p.Results;
if any(strcmp('digits',p.UsingDefaults))
    r.digits=max(floor(log(r.num)/log(r.base)),0)+1;
end

binvec=mod(floor(r.num./r.base.^(0:(r.digits-1))),2);

end

