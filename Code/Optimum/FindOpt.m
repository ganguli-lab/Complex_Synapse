function [ Wp,Wm ] = FindOpt( t,n,varargin )
%[Wp,Wm]=FINDOPT(t,n) Find synapse model that maximises SNR(t)
%   t = time value
%   n = #states
%   Wp = potentiation transition rates
%   Wm = depression transition rates

InitRand=true;
Triangular=false;
varargin=assignApplicable(varargin);

if InitRand
    if Triangular
        W=RandTrans(n);
        w=BinaryWeights(n);
        [Wp,Wm]=TriangleDcmp(W,0.5,w,t);
    else
        Wp=RandTrans(n);
        Wm=RandTrans(n);
    end
else
    [Wp,Wm]=MakeSMS(ones(1,n-1));
%     [Wp,Wm]=DiffJump(n);
end

try
    [Wp,Wm] = ModelOpt( Wp,Wm,t,varargin{:});
catch ME
    return;
end

end

