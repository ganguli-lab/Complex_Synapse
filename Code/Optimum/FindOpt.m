function [ Wp,Wm ] = FindOpt( t,n,varargin )
%[Wp,Wm]=FINDOPT(t,n) Find synapse model that maximises SNR(t)
%   t = time value
%   n = #states
%   Wp = potentiation transition rates
%   Wm = depression transition rates

InitRand=true;
varargin=assignApplicable(varargin);

if InitRand
    Wp=RandTrans(n);
    Wm=RandTrans(n);
else
    [Wp,Wm]=MakeSMS(ones(1,n-1));
%     [Wp,Wm]=DiffJump(n);
end

[ Wp,Wm ] = ModelOpt( Wp,Wm,t,varargin{:});

end

