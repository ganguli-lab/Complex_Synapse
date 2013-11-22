function [ Wp,Wm,snr ] = FindOpt( t,n,reps,varargin )
%[Wp,Wm]=FINDOPT(t,n,reps) Find synapse model that maximises SNR(t)
%   t    = time value
%   n    = #states
%   reps = number of attempts we max over
%   Wp = potentiation transition rates
%   Wm = depression transition rates

existsAndDefault('reps',1);

if reps==1

    InitRand=true;
    Triangular=false;
    varargin=assignApplicable(varargin);

    w=BinaryWeights(n);

    if InitRand
        if Triangular
            W=RandTrans(n);
            [Wp,Wm]=TriangleDcmp(W,0.5,w,t);
            Wp=Stochastify(Wp+eye(n));
            Wm=Stochastify(Wm+eye(n));
        else
            Wp=RandTrans(n);
            Wm=RandTrans(n);
        end
    else
    %    [Wp,Wm]=MakeSMS(ones(1,n-1));
         [Wp,Wm]=DiffJump(n);
    end

    try
        [Wp,Wm,snr] = ModelOpt( Wp,Wm,t,varargin{:});
    %     [Wp,Wm,snr] = GradientModelOpt( Wp,Wm,t,varargin{:});
    catch ME
        snr=SNRcurve(t,Wp,Wm,0.5,w,'UseExpM',true);
        disp(ME.message);
        disp(['In function: ' ME.stack(1).name ' at line ' int2str(ME.stack(1).line) ' of ' ME.stack(1).file]);
        return;
    end

else
    
    DispReps=false;
    varargin=assignApplicable(varargin);
    Wp=[];
    Wm=[];
    snr=0;
    for i=1:reps
        [ Wpt,Wmt,St ] = FindOpt( t,n,1,varargin{:} );
        if St>snr
            Wp=Wpt;
            Wm=Wmt;
            snr=St;
        end
        if DispReps
            disp([int2str(i) '/' int2str(reps)]);
        end
    end
        
end

end