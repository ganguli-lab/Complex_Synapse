function [ qv,A ] = FindOptChainHomLA( s,n,varargin )
%[qv,A]=FINDOPTCHAINHOMLA(s,n,reps) Find synapse model that maximises SNR(t)
%   s    = inverse time value
%   n    = #states
%   reps = number of attempts we max over
%   qv = nearest neighbour transition rates
%   A  = Laplace Transf value

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='FindOptChainHomLA';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addOptional('reps',1,@(x)validateattributes(x,{'numeric'},{'scalar'},'FindOptChainHomLA','reps',3));
    p.addParameter('InitRand',true,@(x) validateattributes(x,{'logical'},{'scalar'},'FindOptChainHomLA','InitRand'));
    p.addParameter('InitHomZero',false,@(x) validateattributes(x,{'logical'},{'scalar'},'FindOptChainHomLA','InitHomZero'));
    p.addParameter('DispReps',false,@(x) validateattributes(x,{'logical'},{'scalar'},'FindOptChainHomLA','InitRand'));
end
p.parse(varargin{:});
r=p.Results;

if r.reps==1

    if r.InitRand
        qv=rand(1,2*n-2);
    else
    %    [Wp,Wm]=MakeSMS(ones(1,n-1));
         qv=ones(1,2*n-2);
    end
    
    if ~r.InitHomZero
        qv=qv+rand(1,2*n-2);
    end

    try
    %     [Wp,Wm] = ModelOpt( Wp,Wm,t,varargin{:});
        [ qv,A ]=ModelOptChainHomLA( qv,s,p.Unmatched);
    catch ME
        A=-OptFunChainHomLA(qv,s,0.5);
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
        [ qvt,At ] = FindOptChainHomLA( s,n,1, p.Unmatched,'InitRand',r.InitRand );
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
        [ qvt,At ] = FindOptChainHomLA( s,n,1, p.Unmatched,'InitRand',r.InitRand,'InitHomZero',true );
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

