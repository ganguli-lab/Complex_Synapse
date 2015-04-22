function [ qv,A ] = FindOptChainHomL( s,n,varargin )
%[Wp,Wm]=FINDOPTCHAINHOML(t,n,reps) Find synapse model that maximises SNR(t)
%   t    = time value
%   n    = #states
%   reps = number of attempts we max over
%   Wp = potentiation transition rates
%   Wm = depression transition rates

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='FindOptChainL';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addOptional('reps',1,@(x)validateattributes(x,{'numeric'},{'scalar'},'FindOptChainL','reps',3));
    p.addParameter('InitRand',true,@(x) validateattributes(x,{'logical'},{'scalar'},'FindOptChainL','InitRand'));
    p.addParameter('DispReps',false,@(x) validateattributes(x,{'logical'},{'scalar'},'FindOptChainL','InitRand'));
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

    try
    %     [Wp,Wm] = ModelOpt( Wp,Wm,t,varargin{:});
        [ qv,A ]=ModelOptChainHomL( qv,s,p.Unmatched);
    catch ME
        A=-OptFunChainHomL(qv,s);
        disp(ME.message);
        disp(['In function: ' ME.stack(1).name ' at line ' int2str(ME.stack(1).line) ' of ' ME.stack(1).file]);
        return;
    end

else

    warning('off','MATLAB:nearlySingularMatrix');
    qv=[];
    A=0;
    for i=1:r.reps
        [ qvt,At ] = FindOptChainHomL( s,n,1, p.Unmatched,'InitRand',r.InitRand );
        if At>A
            qv=qvt;
            A=At;
        end
        if r.DispReps
            disp([int2str(i) '/' int2str(r.reps)]);
        end
    end
    warning('on','MATLAB:nearlySingularMatrix');
end

end

