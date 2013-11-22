function [ qv,A ] = FindOptChainL( s,n,reps,varargin )
%[Wp,Wm]=FINDOPT(t,n,reps) Find synapse model that maximises SNR(t)
%   t    = time value
%   n    = #states
%   reps = number of attempts we max over
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
        A=-OptFunChainL(qv,s);
        disp(ME.message);
        disp(['In function: ' ME.stack(1).name ' at line ' int2str(ME.stack(1).line) ' of ' ME.stack(1).file]);
        return;
    end

else

    DispReps=false;
    varargin=assignApplicable(varargin);
    qv=[];
    A=0;
    for i=1:reps
        [ qvt,At ] = FindOptChainL( s,n,1,varargin{:} );
        if At>A
            qv=qvt;
            A=At;
        end
        if DispReps
            disp([int2str(i) '/' int2str(reps)]);
        end
    end
    
end

end

