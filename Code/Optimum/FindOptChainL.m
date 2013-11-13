function [ qv,A ] = FindOptChainL( s,n,reps,varargin )
%[Wp,Wm]=FINDOPT(t,n) Find synapse model that maximises SNR(t)
%   t = time value
%   n = #states
%   Wp = potentiation transition rates
%   Wm = depression transition rates

existsAndDefault('reps',1);

if reps==1

    InitRand=true;
    varargin=assignApplicable(varargin);

    if InitRand
        qv=rand(1,n-1);
    else
    %    [Wp,Wm]=MakeSMS(ones(1,n-1));
         qv=ones(1,n-1);
    end

    try
    %     [Wp,Wm] = ModelOpt( Wp,Wm,t,varargin{:});
        [ qv,A ]=ModelOptChainL( qv,s,varargin{:});
    catch ME
        A=OptFunChainL(qv,s);
        return;
    end

else

    qv=[];
    A=0;
    for i=1:reps
        [ qvt,At ] = FindOptChainL( s,n,1,varargin{:} );
        if At>A
            qv=qvt;
            A=At;
        end
    end
    
end

end

