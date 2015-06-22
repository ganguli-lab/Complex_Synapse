function [ Wp,Wm,A ] = FindOptDoubleL( sm,sc,Ac,n,varargin )
%[Wp,Wm,A]=FINDOPTDOUBLEL(sm,n,reps) Find synapse model that maximises
%A(sm) subject to A(sc)=Ac
%   sm   = inverse time value
%   n    = #states
%   reps = number of attempts we max over
%   Wp = potentiation transition rates
%   Wm = depression transition rates
%   A  = Laplace Transf value

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='FindOptL';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addOptional('reps',1,@(x)validateattributes(x,{'numeric'},{'scalar'},'FindOptL','reps',3));
    p.addParameter('InitRand',true,@(x) validateattributes(x,{'logical'},{'scalar'},'FindOptL','InitRand'));
%     p.addParameter('Triangular',false,@(x) validateattributes(x,{'logical'},{'scalar'},'FindOpt','InitRand'));
    p.addParameter('DispReps',false,@(x) validateattributes(x,{'logical'},{'scalar'},'FindOptL','InitRand'));
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
    %     end
    else
    %    [Wp,Wm]=MakeSMS(ones(1,n-1));
         [Wp,Wm]=DiffJump(n);
    end

    try
    %     [Wp,Wm] = ModelOpt( Wp,Wm,t,varargin{:});
        [Wp,Wm,A] = ModelOptDoubleL( Wp,Wm,sm,sc,Ac,p.Unmatched);
    catch ME
        A=SNRlaplace(sm,Wp,Wm,0.5,BinaryWeights(n));
        disp(ME.message);
        disp(['In function: ' ME.stack(1).name ' at line ' int2str(ME.stack(1).line) ' of ' ME.stack(1).file]);
        return;
    end
    
else
    
    warning('off','MATLAB:nearlySingularMatrix');
    Wp=[];
    Wm=[];
    A=0;
    for i=1:r.reps
        if r.DispReps
            DispCounter(i,r.reps,'rep: ');
        end
%         [ Wpt,Wmt,At ] = FindOptL( sm,n,1, p.Unmatched,'InitRand',r.InitRand,'Triangular',r.Triangular );
        [ Wpt,Wmt,At ] = FindOptL( sm,n,1, p.Unmatched,'InitRand',r.InitRand );
        if At>A
            Wp=Wpt;
            Wm=Wmt;
            A=At;
        end
    end
        if r.DispReps
            DispCounter(r.reps+1,r.reps,'rep: ');
        end
   
    warning('on','MATLAB:nearlySingularMatrix');
end

end

