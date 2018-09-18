function [ tf,ix ] = istransient( W,varargin )
%[tf,ix]=ISTRANSIENT(W,thresh) Is Markov chain transient?
%   tf = are there any transient states?(true/false
%   ix = indices of transient states
%   W = transition rate matrix
%   thresh = threshold for transience, sum_i(i~=j) W_ij < thresh, (default: 1e-5)

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='istransient';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addRequired('W',@(x)validateattributes(x,{'numeric'},{'2d','square'},'istransient','W',1))
    p.addOptional('thresh',1,@(x)validateattributes(x,{'numeric'},{'scalar','positive'},'istransient','thresh',2));
    p.addParameter('UseP',true,@(x) validateattributes(x,{'logical'},{'scalar'},'istransient','UseP'));
end
p.parse(W,varargin{:});
r=p.Results;



if r.UseP
    pr=EqProb(r.W);
    tf= pr <r.thresh;
else
    tf = sum(r.W,1)-diag(r.W)' < thresh;
end

ix=find(tf);
tf=any(tf);


end

