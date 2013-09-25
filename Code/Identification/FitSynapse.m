function [ fitmodel,loglike,exitflag,output ] = FitSynapse( simobj,guessmodel,options )
%FITSYNAPSE Summary of this function goes here
%   Detailed explanation goes here

defoptions=struct('MaxIter',1000,'TolFun',-1e-6,'TolX',1e-6,'TolFunChange',1e-6,...
    'Algorithm','BW','Weighter','RJ','ExtraOpts',{},...
    'Display','off','OutputFcn',[],'PlotFcn',[],...
    'GroundTruth',[]);
if isscalar(simobj)
    defoptions.Weighter='Uni';
end
defoptions=UpdateOptions(defoptions,options);

fitmodel=guessmodel.Sort;
loglike=HMMloglike(fitmodel,simobj);

weighter=str2func([defoptions.Weighter 'weight']);

exitflag=0;
msg=['Exceeded max iterations: ' int2str(defoptions.MaxIter)];

if ~isempty(defoptions.GroundTruth)
    trueloglike=HMMloglike(defoptions.GroundTruth,simobj);
end

for i=1:defoptions.MaxIter
    prevmodel=fitmodel;
    prevloglike=loglike;
    [fitmodel,~,loglike]=weighter(prevmodel,simobj,'Algorithm',defoptions.Algorithm,defoptions.ExtraOpts{:});
    [kli,klm]=prevmodel.KLDivs(fitmodel);
    %
    if ~isempty(defoptions.GroundTruth)
        [truekli,trueklm]=defoptions.GroundTruth.KLDivs(fitmodel);
        if trueloglike-loglike < defoptions.TolFun
            if mean(trueklm) < defoptions.TolX
                if truekli < defoptions.TolX
                    exitflag=1;
                    msg=['Success. trueloglike - loglike < ' num2str(defoptions.TolFun)...
                        ' and KLdiv from true model to fit model < ' num2str(defoptions.TolX)];
                    break;
                else
                    exitflag=-4;
                    msg=['Not enough data? trueloglike - loglike < ' num2str(defoptions.TolFun)...
                        ' and KLdiv from true M to fit M > ' num2str(defoptions.TolX)...
                        ' but not Initial'];
                    break;
                end
            else
                exitflag=-3;
                msg=['Not enough data? trueloglike - loglike < ' num2str(defoptions.TolFun)...
                    ' despite KLdiv from true M to fit M > ' num2str(defoptions.TolX)];
                break;
            end
        end
    end
    %
    if loglike > defoptions.TolFun
        exitflag=1;
        msg=['Success. loglike > ' num2str(defoptions.TolFun)];
        break;
    elseif mean([kli klm]) < defoptions.TolX && abs(loglike-prevloglike) < defoptions.TolFunChange
        exitflag=-2;
        msg=['Reached local maximum. Change in loglike < ' num2str(defoptions.TolFunChange)...
            '. KL-div due to change in model < ' num2str(defoptions.TolX)];
        break;
    end
end


output=struct('message',msg,'algortihm',defoptions.Algorithm,'weighter',defoptions.Weighter,...
    'iterations',i,'KLstepInitial',kli,'KLstepM',klm);

if ~isempty(defoptions.GroundTruth)
    output.groundTruth=defoptions.GroundTruth;
    output.trueloglike=trueloglike;
    output.KLtrueInitial=truekli;
    output.KLtrueM=trueklm;
end

    function opt=UpdateOptions(oldopt,newopt)
        opt=oldopt;
        pr=fields(newopt);
        for ii=1:length(pr)
            if isfield(opt,pr{ii}) && ~isempty(newopt.(pr{ii}))
                opt.(pr{ii})=newopt.(pr{ii});
            end
        end
    end

end

