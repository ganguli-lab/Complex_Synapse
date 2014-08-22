function [ genfun ] = MultiExpGenFunAC( s,data,varargin )
%genfun=MULTIEXPGENFUNAC(s,data) analytic continuation of estimate of
%moment generating function
%   genfun = <exp(-st)> for Re(s)>0
%   s      = column of s values
%   data   = row of ssamples of t values

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='MultiExpGenFunAC';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addParameter('R',100);
    p.addParameter('ds',0.1);
    p.addParameter('eps',1);
end
p.parse(varargin{:});
r=p.Results;


s_cont = (-r.R:r.ds:r.R)*1i+r.eps;

genfun_cont = MultiExpGenFun(s_cont',data);

kernel = 1./(bsxfun(@minus,s_cont,s));

genfun = kernel*genfun_cont*r.ds/(2*pi);


end

