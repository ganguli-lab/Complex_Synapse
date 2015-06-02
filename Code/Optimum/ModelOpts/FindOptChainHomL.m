function [ qv,A ] = FindOptChainHomL( s,n,varargin )
%[qv,A]=FINDOPTCHAINHOML(s,n,reps) Find synapse model that maximises SNR(t)
%   s    = inverse time value
%   n    = #states
%   reps = number of attempts we max over
%   qv = nearest neighbour transition rates
%   A  = Laplace Transf value

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='FindOptChainHomL';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addOptional('reps',1,@(x)validateattributes(x,{'numeric'},{'scalar'},'FindOptChainHomL','reps',3));
    p.addParameter('InitRand',true,@(x) validateattributes(x,{'logical'},{'scalar'},'FindOptChainHomL','InitRand'));
    p.addParameter('InitHomZero',false,@(x) validateattributes(x,{'logical'},{'scalar'},'FindOptChainHomL','InitHomZero'));
    p.addParameter('DispReps',false,@(x) validateattributes(x,{'logical'},{'scalar'},'FindOptChainHomL','InitRand'));
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
        qv(2*n-1:end)=0;
    end

    try
    %     [Wp,Wm] = ModelOpt( Wp,Wm,t,varargin{:});
        [ qv,A ]=ModelOptChainHomL( qv,s,p.Unmatched);
    catch ME
        A=-OptFunChainHomL(qv,s,0.5);
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
        [ qvt,At ] = FindOptChainHomL( s,n,1, p.Unmatched,'InitRand',r.InitRand );
        if At>A
            qv=qvt;
            A=At;
        end
        if r.DispReps
            DispCounter(i,r.reps,'rep: ');
%             disp([int2str(i) '/' int2str(r.reps)]);
        end
    end
            DispCounter(r.reps+1,r.reps,'rep: ');
    DispCounter(2,2,'init: ');
    for i=1:r.reps
        [ qvt,At ] = FindOptChainHomL( s,n,1, p.Unmatched,'InitRand',r.InitRand,'InitHomZero',true );
        if At>A
            qv=qvt;
            A=At;
        end
        if r.DispReps
            DispCounter(i,r.reps,'rep: ');
%             disp([int2str(i) '/' int2str(r.reps)]);
        end
    end
            DispCounter(r.reps+1,r.reps,'rep: ');
    DispCounter(3,2,'init: ');
    warning('on','MATLAB:nearlySingularMatrix');
end

end

