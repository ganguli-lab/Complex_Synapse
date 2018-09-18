function [ S ] = SNRcurveWq( t, W, q ,fp, w, varargin  )
%S=SNRCURVE(T,W,q,fp,w) SNR as function of time
%   T = time values
%   W = forgetting transition rates: f^+W^+ + f^-W^-
%   q = encoding transition rates: W^+ - W^-
%   w = Weights of states (+/-1)

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='SNRcurveWq';
    p.CaseSensitive=true;
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addRequired('t',@(x)validateattributes(x,{'numeric'},{},'SNRcurveWq','t',1));
    p.addRequired('W',@(x)validateattributes(x,{'numeric'},{'2d','square'},'SNRcurveWq','W',2));
    p.addRequired('q',@(x)validateattributes(x,{'numeric'},{'2d','square'},'SNRcurveWq','q',3));
    p.addRequired('fp',@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative','<=',1},'SNRcurveWq','fp',4));
    p.addRequired('w',@(x)validateattributes(x,{'numeric'},{'column'},'SNRcurveWq','w',5));
    p.addParameter('UseExpM',false,@(x)validateattributes(x,{'logical'},{'scalar'},'SNRcurveWq','UseExpM'));
end
p.parse(t,W,q,fp,w,varargin{:});
r=p.Results;

error(CheckSize(r.q,@(x)samesize(r.W,x),'samesize(W)'));%also square matrix of same size
error(CheckSize(r.w,@(x)length(x)==length(r.W),'samesize(W)'));

if r.UseExpM
    S=zeros(size(r.t));
    p=EqProb(r.W,p.Unmatched);
    for i=1:numel(r.t)
        S(i) = p*r.q*expm(r.W*r.t(i))*r.w;
    end
    S=2*r.fp*(1-r.fp)*S;
else
    [ qa,ca ] = SpectrumWq(r.W,r.q,r.fp,r.w,p.Unmatched);
    S = SNRcurveCaQa(r.t,ca,qa);
end

end

