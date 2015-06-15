function [ qv,A ] = FindOptChainHomLC( s,n,varargin )
%[qv,A]=FINDOPTCHAINHOMLC(s,n,reps) Find synapse model that maximises SNR(t)
%   s    = inverse time value
%   n    = #states
%   reps = number of attempts we max over
%   Wp = potentiation transition rates
%   Wm = depression transition rates

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='FindOptChainHomLC';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addOptional('reps',1,@(x)validateattributes(x,{'numeric'},{'scalar'},'FindOptChainHomLC','reps',3));
    p.addParameter('InitRand',true,@(x) validateattributes(x,{'logical'},{'scalar'},'FindOptChainHomLC','InitRand'));
    p.addParameter('InitHomZero',false,@(x) validateattributes(x,{'logical'},{'scalar'},'FindOptChainHomLC','InitHomZero'));
    p.addParameter('DispReps',false,@(x) validateattributes(x,{'logical'},{'scalar'},'FindOptChainHomLC','InitRand'));
end
p.parse(varargin{:});
r=p.Results;

if r.reps==1

    if r.InitRand
        qv=rand(1,4*n-4);
    else
    %    [Wp,Wm]=MakeSMS(ones(1,n-1));
         qv=ones(1,4*n-4);
    end
    
    if r.InitHomZero
        qv(n:3*n-3)=0;
    end

    try
    %     [Wp,Wm] = ModelOpt( Wp,Wm,t,varargin{:});
        [ qv,A ]=ModelOptChainHomLC( qv,s,p.Unmatched);
    catch ME
        A=-OptFunChainHomLC(qv,s,0.5);
        disp(ME.message);
        disp(['In function: ' ME.stack(1).name ' at line ' int2str(ME.stack(1).line) ' of ' ME.stack(1).file]);
        return;
    end

else

    warning('off','MATLAB:nearlySingularMatrix');
    qv=[];
    A=0;
    DispCounter(1,2,'init: ');
    for i=1:r.reps
        [ qvt,At ] = FindOptChainHomLC( s,n,1, p.Unmatched,'InitRand',r.InitRand );
        if At>A
            qv=qvt;
            A=At;
        end
        if r.DispReps
            DispCounter(i,r.reps,'rep: ');
%             disp([int2str(i) '/' int2str(r.reps)]);
        end
    end
        if r.DispReps
            DispCounter(r.reps+1,r.reps,'rep: ');
        end
    DispCounter(2,2,'init: ');
    for i=1:r.reps
        [ qvt,At ] = FindOptChainHomLC( s,n,1, p.Unmatched,'InitRand',r.InitRand,'InitHomZero',true );
        if At>A
            qv=qvt;
            A=At;
        end
        if r.DispReps
            DispCounter(i,r.reps,'rep: ');
%             disp([int2str(i) '/' int2str(r.reps)]);
        end
    end
        if r.DispReps
            DispCounter(r.reps+1,r.reps,'rep: ');
        end
    DispCounter(3,2,'init: ');
    warning('on','MATLAB:nearlySingularMatrix');
end

end

