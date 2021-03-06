function [ Wp,Wm,Q,A ] = FindOptHomL( sm,n,varargin )
%[Wp,Wm,Q,A]=FINDOPTHOML(sm,n,reps) Find synapse model that maximises A(sm)
%   sm   = inverse time value
%   n    = #states
%   reps = number of attempts we max over
%   Wp = potentiation transition rates
%   Wm = depression transition rates
%   Q  = activity independent transition rates
%   A  = Laplace Transf value

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='FindOptHomL';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addOptional('reps',1,@(x)validateattributes(x,{'numeric'},{'scalar'},'FindOptHomL','reps',3));
    p.addParameter('InitRand',true,@(x) validateattributes(x,{'logical'},{'scalar'},'FindOptHomL','InitRand'));
    p.addParameter('InitHomZero',false,@(x) validateattributes(x,{'logical'},{'scalar'},'FindOptHomL','InitHomZero'));
%     p.addParameter('Triangular',false,@(x) validateattributes(x,{'logical'},{'scalar'},'FindOptHomL','InitRand'));
    p.addParameter('DispReps',false,@(x) validateattributes(x,{'logical'},{'scalar'},'FindOptHomL','InitRand'));
end
p.parse(varargin{:});
r=p.Results;

if r.reps==1


    if r.InitRand
    %     if r.Triangular
    %         W=RandTrans(n);
    %         w=BinaryWeights(n);
    %         [Wp,Wm]=TriangleDcmp(W,0.5,w,sm);
    %         Wp=Stochastify(Wp+eye(n));
    %         Wm=Stochastify(Wm+eye(n));
    %     else
            Wp=RandTrans(n);
            Wm=RandTrans(n);
            Q=RandTrans(n);
    %     end
    else
    %    [Wp,Wm]=MakeSMS(ones(1,n-1));
         [Wp,Wm]=DiffJump(n);
         Q=zeros(n);
    end

    if r.InitHomZero
        Q=zeros(n);
    end
    
    
    try
    %     [Wp,Wm] = ModelOpt( Wp,Wm,t,varargin{:});
        [Wp,Wm,Q,A] = ModelOptHomL( Wp,Wm,Q,sm,p.Unmatched);
    catch ME
        A=SNRlaplace(sm,Wp+Q,Wm+Q,0.5,BinaryWeights(n));
        disp(ME.message);
        disp(['In function: ' ME.stack(1).name ' at line ' int2str(ME.stack(1).line) ' of ' ME.stack(1).file]);
        return;
    end
    
else
    warning('off','MATLAB:nearlySingularMatrix');
    
    Wp=[];
    Wm=[];
    Q=[];
    A=0;
    DispCounter(1,2,'init: ');
    for i=1:r.reps
        if r.DispReps
            DispCounter(i,r.reps,'rep: ');
        end
%         [ Wpt,Wmt,At ] = FindOptL( sm,n,1, p.Unmatched,'InitRand',r.InitRand,'Triangular',r.Triangular );
        [ Wpt,Wmt,Qt,At ] = FindOptHomL( sm,n,1, p.Unmatched,'InitRand',r.InitRand );
        if At>A
            Wp=Wpt;
            Wm=Wmt;
            Q=Qt;
            A=At;
        end
    end
        if r.DispReps
            DispCounter(r.reps+1,r.reps,'rep: ');
        end
    DispCounter(2,2,'init: ');
    for i=1:r.reps
        if r.DispReps
            DispCounter(i,r.reps,'rep: ');
        end
%         [ Wpt,Wmt,At ] = FindOptL( sm,n,1, p.Unmatched,'InitRand',r.InitRand,'Triangular',r.Triangular );
        [ Wpt,Wmt,Qt,At ] = FindOptHomL( sm,n,1, p.Unmatched,'InitRand',r.InitRand,'InitHomZero',true );
        if At>A
            Wp=Wpt;
            Wm=Wmt;
            Q=Qt;
            A=At;
        end
    end
        if r.DispReps
            DispCounter(r.reps+1,r.reps,'rep: ');
        end
    DispCounter(3,2,'init: ');
    warning('on','MATLAB:nearlySingularMatrix');
    
end

end

