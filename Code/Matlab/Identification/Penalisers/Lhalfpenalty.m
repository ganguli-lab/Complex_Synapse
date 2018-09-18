function [ newmodelobj,loglike,pstate ] = Lhalfpenalty( modelobj,simobj,varargin )
%[newmodelobj,loglike,pstate]=OFFDIAGL1PENALTY(modelobj,simobj) update of
%estiamted HMM L^1/2 penalty
%   newmodelobj = updated SynapseIdModel
%   loglike     = log likelihood of readouts under old model (prod over simobj)
%   pstate      = posterior prob of HMM being in each state at each time (cell, one element for each simobj)
%   modelobj = SynapseIdModel
%   simobj   = vector of SynapsePlastSeq

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='Lhalfpenalty';
    p.StructExpand=true;
    p.KeepUnmatched=true;
    p.addParameter('Weighter','RJ');
    p.addParameter('Normalise',true);
    p.addParameter('Penalty',1);
%     p.addParameter('Normalise',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
end
p.parse(varargin{:});
Penalty=p.Results.Penalty;

weighterfn=str2func([p.Results.Weighter 'weight']);

[newmodelobj,loglike,pstate] = weighterfn(modelobj,simobj,'Normalise',false,p.Unmatched);

if Penalty>0
    M=newmodelobj.M;
    L1=0;
    for i=1:newmodelobj.NumPlast
        M{i}=M{i}/Penalty^2;
        for j=1:modelobj.NumStates
            gamma=CalcLagrange(M{i}(j,:));
            M{i}(j,:)=PenNorm(M{i}(j,:),gamma);
        end%for j
        M{i}=M{i}*Penalty^2;
        L1=L1+sum(sqrt(M{i}(:)));
    end%for i
    newmodelobj=newmodelobj.setM(M);

    loglike = loglike - Penalty * L1;
end%if Penalty

if p.Results.Normalise
    newmodelobj=newmodelobj.Normalise;
    assert(newmodelobj.isvalid,'newmodelobj is invalid');
end

    function gamma=CalcLagrange(Lrow)
        %find Lagrange multiplier for each row
        %   Lrow=Mrow/beta^2
        interval=[0 100*(2/length(Lrow)*Penalty^2)^(1/4)];
        gamma = fzero( @(x) LagrangeCond(x,Lrow), interval );
    end

    function val=LagrangeCond(gamma,Lrow)
        %correct value of Lagrange multiplier given by root of this function
        %   Lrow=Mrow/beta^2
        %   gamma^2 = 1/lambda
        val = sum( PenNorm(Lrow,gamma) ) - 1 / Penalty^2 ;
    end

    function Lnew=PenNorm(Lrow,gamma)
        %   Lnew=Mnew/beta^2
        %   Lrow=Mrow/beta^2
        Lnew = ( 2*gamma^2 * Lrow + gamma^4 - gamma^3 * sqrt( gamma^2 + 4*Lrow ) ) / 2 ;
    end

end

